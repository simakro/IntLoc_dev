# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

"""
 This module evaluates the probability for other potential integration sites to
 exist in the genome. This is achieved by copying the sequences around the already
 identified integration sites, which are covered by evidence reads, and using them
 as query against the reference. The first hit will be the sequence itself, and
 the second best hit is the one to compare against. Theoretically, the self hit 
 could be prevented by cropping the primary int-site from the reference instead
 of just copying it, but this strategy would fail, in cases where there are one or
 more highly identical sequences to the real integration site and the reads would
 likely be be distributed more or less equally between these different sites during
 mapping. Then, when generating the cropped/masked sequence, all the other potential 
 candidates would be cropped as well, because they co-exist as primary identified 
 int-sites. This would lead to all other misleading candidates to be removed and 
 thus defeat the purpose, because it would especially happen if there are very similar
 potential integration sites.
 Intloc offer a number of metrics: percentage of integration site covered by reads
 (sum qcovhsp/qcovus); edit distance/identity; read coverage (RE); 
 bitscore/evalues (cumulative or similar); alternative integration sites 
"""

import os
import sys
import json
from collections import defaultdict
import subprocess
import multiprocessing as mp
import pandas as pd
from copy import copy

from il_logging import Ilogger
from il_chk_conv import get_blast_stats, check_input_file, report_time, get_reads_by_header
from il_classes import LocalAlnm, ContAlnm, ReadToRefAl
from il_asm_db_tools import gen_blastdb, cleanup_blastdb
from il_subprocesses import minimap2_sp, caplog_sp_error
import il_map_figs


ilog = Ilogger()
ilog.module = __name__

def load_alignments(outdir):
    """Load alignments from json file containing candidate reads mapped to ref 
    exported from PotIntSite objects"""
    ilog.vlprint("Loading alignments", 3)
    intsite_als = os.path.join(outdir, "IntSitesAlignments.json")
    with open(intsite_als, "r") as alj:
        alignments = json.load(alj)
    return alignments


class IntSite(ReadToRefAl):
    def __init__(self, name, primary=True):
        self.name = name
        self.old_name = False
        self.assoc_primary_site = ""
        self.primary_site = primary
        self.chrom = "_".join(name.split("_")[:-1]) if self.primary_site else ""
        self.pos = int(name.split("_")[-1]) if self.primary_site else ""
        self.alignments = []
        self.non_int_alignments = []
        self.spanned_region = ()
        self.spanned_length = int()
        self.int_pos_in_span_reg = int()
        self.spanned_seq = ""
        self.similar_regions = []
        self.contig_sim_regs = []
        self.eval_alignments = {}
        self.eval_cum_bitscore = int()


    def analyze_similar_regions(self):
        """Select most similar region"""
        if len(self.similar_regions) > 0:
            self.most_similar = max(self.similar_regions, key=lambda x: x.bitscore)
            bitscores = [hsp.bitscore for hsp in self.similar_regions]
            # below attributes may become useful at a later point
            self.sum_simil_bitscores = sum(bitscores)
            self.avg_simil_bitscores = self.sum_simil_bitscores / len(bitscores)
            return True
        else:
            return False

    def get_spanned_region(self):
        """Determine genomic region spanned by mapped candidate reads"""
        corner_vals = [val for list in self.alignments for val in list]
        corner_vals = sorted([val for list in corner_vals for val in list])
        max_span = (corner_vals[0], corner_vals[-1])
        self.spanned_region = max_span
        self.spanned_length = max_span[1] - max_span[0] + 1
        self.int_pos_in_span_reg = self.pos - max_span[0]
    
    # def get_spanned_seq(self, chrom):
    #     """Extract spanned sequence from reference"""
    #     self.spanned_seq = chrom[self.spanned_region[0]:self.spanned_region[1]+1]
    
    def generate_contiguous_sim_regs(self):
        al_dict = defaultdict(list)
        for loc_al in self.similar_regions:
            al_dict[loc_al.chrom].append(loc_al)
        self.contig_sim_regs = [
            ContAlnm(
                self.spanned_length,
                chrom_key,
                al_dict[chrom_key]
                ) for chrom_key in al_dict
            ]
        for contal in self.contig_sim_regs:
            contal.refine_self()
    
    def analyze_sec_site(self):
        """Analyze secondary site for cumlative raw- and bitscore"""
        locals = defaultdict(list)
        with open(f"eval_{self.assoc_primary_site}_second", "r") as res:
            for line in res:
                if line.startswith("RNAME"):
                    pass
                else:
                    ls = line.split(',')
                    new_local = LocalAlnm(*ls)
                    locals[new_local.read].append(new_local)
        contals = {}
        for read in locals:
            readlen = locals[read][0].read_len
            match_subject = locals[read][0].chrom
            new_contal = ContAlnm(readlen, match_subject, locals[read])
            new_contal.refine_self()
            contals[read] = new_contal
        self.eval_alignments = contals
        self.eval_cum_bitscore = sum(
            [sum(contal.bitscores) for contal in self.eval_alignments.values()]
            )
        return self.eval_cum_bitscore


def read_out_spanned_regions(alignments, primary_int_site_objs):
    """Generate primary IntSite objects and determine spanned regions"""
    ilog.vlprint("Generating primary IntSite objects", 3)
    intsites = {}
    pisos = {site.name:site for site in primary_int_site_objs}
    pisold = {s.old_name:s for s in primary_int_site_objs if s.old_name}
    for intsite in alignments:
        site = IntSite(intsite)
        site.assoc_primary_site = site.name
        for alnm in alignments[intsite]:
            chrom_coord = sorted(
                [sub_al[1] for sub_al in alnm['al_coord']],
                key=lambda x: x[0]
                )
            site.alignments.append(chrom_coord)
        intsites[intsite] = site
    ilog.vlprint("Determining spanned regions", 1)
    for site in intsites:
        site_name = intsites[site].assoc_primary_site
        intsites[site].get_spanned_region()
        try:
            pisos[site_name].spanned_region["minimap2"] = intsites[site].spanned_region
        except KeyError:
            pisold[site_name].spanned_region["minimap2"] = intsites[site].spanned_region
            ilog.vlprint(
                f"Spanned-region was sorted to {pisold[site_name].name} using "
                f"old name {pisold[site_name].old_name}", 3
                )
    return intsites


def crop_spans(req_chroms, cname, spanned_seqs, curr_chrom_lines):
    cchrom_seq = "".join(curr_chrom_lines)
    for stup in req_chroms[cname]:
        start, end = stup[1][0], stup[1][1] + 1
        start = start if start > 0 else 0
        end = end if end <= len(cchrom_seq)+1 else len(cchrom_seq)+1
        spanned_seqs[stup[0]] = cchrom_seq[start:end]


def get_spanned_seqs(intsite_objs, ref_genome):
    """Load required reference chromosomes and extract spanned regions"""
    ilog.vlprint(f"Extracting spanned region", 1)
    site_objs = intsite_objs.values()
    req_chroms = defaultdict(list)

    for site in site_objs:
        req_chroms[site.chrom].append(
            (site.assoc_primary_site, site.spanned_region)
            )

    with open(ref_genome, "r") as ref:
        contigs = {}
        spanned_seqs = {}
        req_chrom = False
        curr_chrom_lines = []
        cname = ""
        for line in ref:
            if line.startswith(">"):
                if len(curr_chrom_lines) > 0:
                    crop_spans(req_chroms, cname, spanned_seqs, curr_chrom_lines)
                    req_chrom = False
                    curr_chrom_lines, cname = [], ""
                cname = line.split(" ")[0].split(">")[1]
                contigs[cname] = 0
                req_chrom = True if cname in req_chroms else False    
            else:
                if req_chrom:
                    curr_chrom_lines.append(line.strip())
        if len(curr_chrom_lines) > 0:
            crop_spans(req_chroms, cname, spanned_seqs, curr_chrom_lines)
    for site in intsite_objs:
        intsite_objs[site].spanned_seq = spanned_seqs[site]

    return intsite_objs


def search_similar_refgen(intsite_objs, ref_genome, args):
    """Search in ref genome for sequences similar to integration site region"""
    # using minimap2 here may speed up for large number of ints and large ref
    ilog.vlprint(
        "Searching reference for sequences similar to integration sites", 1
        )
    temp = os.path.join(args.outdir, "eval_int_sites_spanseq.que")
    out = os.path.join(args.outdir, "eval_int_sites_sim_search.bre")

    if os.path.exists(temp):
        os.remove(temp)
    for int_site in intsite_objs:
        span_seq = intsite_objs[int_site].spanned_seq
        with open(temp, "a") as tmp:
            tmp.write(f">{int_site}_spanned_context\n")
            tmp.write(span_seq + "\n")

    subprocess.run([
                    args.blastn, 
                    "-db", ref_genome,
                    "-query", temp,
                    "-out", out,
                    "-outfmt", "10 qaccver qlen saccver slen length qstart qend"
                    " sstart send pident mismatch gapopen evalue bitscore"
                    " qcovhsp qcovus score gaps nident positive sstrand",
                    "-num_threads", args.cores,
                    "-perc_identity", "90",
                    "-qcov_hsp_perc", "2",
                    "-culling_limit", "5",
                    ],
                    check=True, capture_output=True
                    )
    return out


def filt_similar_seqs(intsite_objs, simil_search_res):
    """Exclude results that lie within the primary intsite region itself, and
        and select longest contiguous alignment"""
    ilog.vlprint("Filter redundant hits within primary region", 4)
    with open(simil_search_res, 'r') as res:
        for line in res:
            ls = line.split(',')
            new_aln = LocalAlnm(*ls)
            qchrom1, qchrom2, intloc = new_aln.read.split('_')[:-2]
            qchrom = '_'.join([qchrom1, qchrom2])
            new_aln.query_chrom = qchrom
            refsite = '_'.join([qchrom, intloc])
            refsite_obj = intsite_objs[refsite]

            if new_aln.chrom == new_aln.query_chrom:
                if new_aln.chrom_start >= refsite_obj.spanned_region[0] + 1:
                    if new_aln.chrom_end <= refsite_obj.spanned_region[1] + 1:
                        pass
            else:
                refsite_obj.similar_regions.append(new_aln)
        for key in intsite_objs:
            iso = intsite_objs[key]
            iso.generate_contiguous_sim_regs()
    return intsite_objs


def check_similar_regions(intsite_objs):
    """Analyze similar regions"""
    ilog.vlprint("Analyzing similar regions", 4)
    for sitename in intsite_objs:
        intsite = intsite_objs[sitename]
        intsite.best_alt_reg = intsite.best_alignment(intsite.contig_sim_regs)
        intsite.analyze_similar_regions()


def create_alt_reg_objs(primary_intsite_objs, ref_genome):
    """Creating alternative IntSite objects where similarity was found"""
    ilog.vlprint("Creating secondary IntSite objects", 4)
    secondary_intsite_objs = {}
    for key in primary_intsite_objs:
        site = primary_intsite_objs[key]
        if site.best_alt_reg:
            sec_name = site.name + "-followup-alt-site"
            sec_site = IntSite(sec_name, primary=False)
            sec_site.assoc_primary_site = site.name
            sec_site.chrom = site.best_alt_reg.chromosome
            al_coord_list = site.best_alt_reg.al_coord
            chrom_coords = [tup[1] for tup in al_coord_list]
            choord_vals = [val for tup in chrom_coords for val in tup]
            sim_range = (min(choord_vals), max(choord_vals))
            # the range of the secondary site similar to the primary site is now 
            # generously extended to both sides by the entire range of the 
            # spanned_seq of the primary site; it is OK if start becomes < 0 or  
            # end > len(chrom) in the process, because this is handled
            # by get_spanned_seqs() internally
            ext_range = (
                sim_range[0]-site.spanned_length, sim_range[0]+site.spanned_length
                )
            sec_site.spanned_region = ext_range
            sec_site.spanned_length = ext_range[1] - ext_range[0] + 1
            # using the primary site name here as key for the secondary site,  
            # so they can later be retrieved from two dicts using the same key
            secondary_intsite_objs[sec_site.assoc_primary_site] = sec_site 
    secondary_intsite_objs = get_spanned_seqs(secondary_intsite_objs, ref_genome)
    
    return secondary_intsite_objs


def blast_eval_secondary(prim_site, blast_db, support_reads, args):
    """Evaluate alternative int-sites by blasting candidate reads against """
    """second best potential location."""
    ilog.vlprint(
        f'Starting BLAST-search to evaluate alternative int-site candidate '
        f'second best to {prim_site}', 3
        )
    os.chdir(args.outdir)
    eval_hits = os.path.join(args.outdir, f"eval_{prim_site}_second")

    try:
        subprocess.run(
                        [args.blastn,
                        "-db", blast_db,
                        "-query", support_reads,
                        "-out", eval_hits,
                        "-outfmt",
                        "10 qaccver qlen saccver slen length qstart qend sstart"
                        " send pident mismatch gapopen evalue bitscore qcovhsp"
                        " qcovus score gaps nident positive sstrand",
                        "-num_threads", "1",
                        "-culling_limit", "1",
                        "-evalue", "0.0001",
                        "-perc_identity", "80"],
                        check=True, capture_output=True
                        )
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, f"BLAST search for evaluation of integration-site with second "
            f"best alternative site for {prim_site} failed", exit=False
            )

    return args.outdir, eval_hits, blast_db


def blast_cand_reads_vs_sec_site_reg(secondary_intsite_objs, args):
    """ 
    1. write sec_span_seqs to seperate fasta files and make-blast-DBs
    2. retrieve cand_reads associated with primary sites
    3. blast cand_reads against sec-site-db
    4. calculate sum bitscores
    5. compare to sum bitscores primary site
    6. Calculate SSBR metric
     """
    ilog.vlprint(f'Starting BLAST-searches for evaluation of int-sites', 1)

    db_jobs = []
    blj_info = {}
    bl_jobs = []
    error_log = open("sec_site_blastdb_err.log", "a")
    for site in secondary_intsite_objs:
        sec_site = secondary_intsite_objs[site]
        db_path = os.path.join(args.outdir, "secsitedb_" + site + ".fa")
        with open(db_path, "w") as db:
            db.write(
                f">{sec_site.name} {sec_site.chrom}~{sec_site.spanned_region[0]}"
                f"-{sec_site.spanned_region[1]}\n"
                )
            db.write(sec_site.spanned_seq)
        cif = check_input_file(db_path, "extensive", interact=False, silent=True)
        corrupt = cif[0]
        if not corrupt:
            blj_info[site] = db_path
            kw_args = {"msg": "secondary int_site candidate"}
            dbp = mp.Process(target=gen_blastdb, args=(db_path, args), kwargs=kw_args)
            db_jobs.append(dbp)
            dbp.start()
        else:
            error_log.write(site + "\n")
    error_log.close()

    for proc in db_jobs:
        proc.join()

    for site in blj_info:
        supp_reads = site + "_supporting_reads.fasta"
        blp = mp.Process(
                        target=blast_eval_secondary, 
                        args=[
                              site,
                              blj_info[site],
                              supp_reads,
                              args
                              ]
                        )     
        bl_jobs.append(blp)
        blp.start()
    
    for proc in bl_jobs:
        proc.join()
    
    # get blast statistics with first site in blj_info
    if len(list(blj_info.keys())) > 0:
        test_sec_site = list(blj_info.keys())[0]
        db_path = blj_info[test_sec_site]
        stats_out = os.path.join(args.outdir, f"stats_{test_sec_site}_out")
        blast_params = get_blast_stats(supp_reads, stats_out, db_path, args)
    else:
        blast_params = "not available"
    # delete all blast-dbs in outdir-folder
    cleanup_blastdb(args.outdir)
    return blast_params


def chk_cr_vs_ssr_success(primary_int_site_objs):
    "Checking successful creation of sec_site blast dbs"
    ilog.vlprint("Checking successful creation of sec_site blast dbs", 1)
    pisos = {piso.name: piso for piso in primary_int_site_objs}
    with open("sec_site_blastdb_err.log", "r") as errl:
        for line in errl:
            ls = line.strip()
            try:
                pisos[ls].ssbr = "n.d."
            except:
                ilog.vlprint(f"INFO: {ls} not in pisos dict",2)


def make_mm_paf_df(paf, ava=False, quick_mm2=False):
    """Generate pandas dataframe from paf alignment file"""
    ilog.vlprint("Generating pandas dataframe from paf alignment file", 5)
    cols = [
            "qname",
            "qlen",
            "qstart",
            "qend",
            "strand",
            "tname",
            "tlen",
            "tstart",
            "tend",
            "matches",
            "map_len",
            "MAPQ",
            "al_type",
            "minimizer",
            "chain_score_prim",
            "chain_score_sec",
            "seq_div",
            "len_rep_seeds"
            ]

    try:
        data = pd.read_csv(paf, sep="\t", header=None)
        data.columns = cols
        ilog.vlprint(f"loaded {paf}", 10)
    except pd.errors.EmptyDataError as e:
        ilog.vlprint(f"{e} no results in {paf}", 10)
        data = pd.DataFrame([], columns=cols)
    except ValueError as ve:
        ilog.vlprint(f"{ve}", 10)
        ilog.vlprint(f"INFO: Excepting ValueError: {ve}", 10)
        if "Length mismatch" in str(ve):
            ilog.vlprint("INFO: Discarding column chain_score_sec", 3)
            cols.remove("chain_score_sec")
            data = pd.read_csv(paf, sep="\t", header=None)
            data.columns = cols
            ilog.vlprint(f"loaded {paf}", 10)
        else:
            ilog.vlprint(
                "ERROR: Col count in df not matching exp. length of 17-18", 3
                )
            sys.exit()

    return data, cols


def filt_sort_paf_df(mm_df, cols, mm2_mapq=False, no_culling=False):
    """Filter alignment info in pandas dataframe"""
    ilog.vlprint("Filtering alignment info in dataframe", 5)
    filt_mq = mm_df[mm_df["al_type"] == "tp:A:P"]
    mm2_mapq = mm2_mapq if mm2_mapq else 20
    filt_mq = filt_mq[filt_mq["MAPQ"] > mm2_mapq]
    sorted_frame = filt_mq.sort_values(by=['tname', 'tstart'])
    ll = [list(row) for row in sorted_frame[:].to_numpy()]
    clz = [list(zip(cols, alnm)) for alnm in ll]
    mm_alnm = [{k:v for k,v in alnm} for alnm in clz]
    if not no_culling:
        mmal_cull = defaultdict(list)
        for al in mm_alnm:
            mmal_cull[al["qname"]].append(al)
        for read in mmal_cull:
            mmal_cull[read] = sorted(
                mmal_cull[read], key=lambda x: (x["MAPQ"], x["map_len"], x["matches"]),
                    )
            mm_alnm = [v[-1] for v in mmal_cull.values()]
    return mm_alnm


def get_prox_hom_int_alnms(prim_int_site_objs):
    prox_hom = [iso for iso in prim_int_site_objs if len(iso.prox_homolog) > 0]
    all_ph_int_alnms = []
    for iso in prox_hom:
        all_ph_int_alnms.extend(iso.supp_als)
    return all_ph_int_alnms, prox_hom


def mmap_all_reads_to_spanseqs(read_file, prim_int_site_objs, args): # site_objs,
    """Mapping subroutine for mmap_total_cov_primsite"""
    primsites = "eval_int_sites_spanseq.que"
    mm2_out = "eval_total_cov_primsites.paf"
    ilog.vlprint("coverage check for evaluation: mapping reads", 5)
    try:
        # subprocess.run(
        #         [
        #         args.minimap2,
        #         "--secondary=no", # available in minimap2 2.3 and higher
        #         # "-N", "0",#can reduce mapping accuracy (def=5) use --secondary
        #         "-x", f"map-{args.seq_tec}",
        #         primsites,
        #         read_file,
        #         "-o", mm2_out  # available in minimap2 2.16 and higher
        #         ],
        #     check=True, capture_output=True
        #     )
        minimap2_sp(primsites, read_file, mm2_out, args, secondary="no")
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, f"Minimap2 for coverage at int-sites failed for {read_file}",
            exit=False
            )
    
    mm_df, cols = make_mm_paf_df(mm2_out)
    mm_alnm = filt_sort_paf_df(mm_df, cols, mm2_mapq=5, no_culling=True)
    piso_name = {pis.name:pis for pis in prim_int_site_objs}
    piso_chrom :dict= defaultdict(list)
    for site in prim_int_site_objs:
        piso_chrom[site.chrom].append(site)
    # additional filter for maplen to exclude aligmnets of short segments
    mm_alnm = [alnm for alnm in mm_alnm if alnm["map_len"]/alnm["qlen"]>0.9]
    for piso in prim_int_site_objs:
        piso.int_als.update({al.READ: al for al in piso.supp_als})
    for al in mm_alnm:
        tsite = "_".join(al["tname"].split("_")[:-2])
        chr_offset = piso_name[tsite].spanned_region["minimap2"][0]
        al["chr"] = "_".join(tsite.split("_")[:-1])
        al["chr_start"] = al["tstart"] + chr_offset
        al["chr_end"] = al["tend"] + chr_offset
        for piso in piso_chrom[al["chr"]]:
            if  al["chr_start"] < piso.approx_loc < al["chr_end"]:
                if al["qname"] not in piso.int_als:
                    piso.pot_non_int_als[al["qname"]] = al
    # # include int_als supporting prox_hom int-site into pot_non_int_als of site on 
    # # other homolog
    # all_ph_intreads, prox_hom = get_prox_hom_int_alnms(prim_int_site_objs)
    # for piso in prox_hom:
    #     matches = [r for r in all_ph_intreads if r.CHR==piso.chrom]
    #     matches = [r for r in matches if r.READ not in piso.int_als]
    #     matches = [
    #         r for r in matches if r.chr_start_end[0]<piso.approx_loc<r.chr_start_end[1]
    #         ]
    #     matches = {r.READ: r for r in matches}
    #     piso.pot_non_int_als.update(matches)
    
    return piso_name


def mmap_total_cov_primsite(all_reads, prim_int_site_objs, args): # site_objs,
    """Determine sequencing coverage at integration sites using minimap2"""
    ilog.vlprint("Determinining sequencing coverage at integration sites", 1)

    if os.path.isdir(all_reads):
        # this takes effect if the args.reads_dir==True
        for file in os.scandir(all_reads):
            piso_dct = mmap_all_reads_to_spanseqs(
                                                    file,
                                                    prim_int_site_objs,
                                                    args
                                                    )
            # for site in als_by_site:
            #     piso_dct[site].pot_non_int_als.update(als_by_site[site]["non_int_reads"])
        for site in piso_dct:
            piso_dct[site].calc_total_cov(mode="potential")
            m1 = len(piso_dct[site].int_als)
            # for read-dir final non_int_read num has to be retrieved from piso
            m2 = len(piso_dct[site].pot_non_int_als)
            m3 = piso_dct[site].total_cov
            ilog.vlprint(
                f"{site} int_als: {m1} pot_non_int_als: {m2} total_cov: {m3}", 2
                )
    else:
        piso_dct = mmap_all_reads_to_spanseqs(
                                                all_reads,
                                                prim_int_site_objs,
                                                args
                                                )
        for site in piso_dct:
            piso_dct[site].calc_total_cov(mode="potential")
            m1 = len(piso_dct[site].int_als)
            m2 = len(piso_dct[site].pot_non_int_als)
            m3 = piso_dct[site].total_cov
            ilog.vlprint(
                f"{site} int_als: {m1} pot_non_int_als: {m2} total_cov: {m3}", 2
                )
            # piso_dct[site].pot_non_int_als = als_by_site[site]["non_int_reads"]


def confirm_non_int_reads(prim_int_site_objs, all_reads, args):
    """Confirm best alignment for non-int-reads"""
    ilog.vlprint("Confirming best alignment for non-int-reads", 2)
    read2loc = defaultdict(dict)
    read2site = defaultdict(list)
    for piso in prim_int_site_objs:
        for read in piso.pot_non_int_als:
            read2site[read].append(piso.name)
            read2loc[read][piso.chrom] = piso.approx_loc

    non_int_file = "non_int_reads.fasta"
    header = list(read2loc.keys())
    # get all genuine non-int-reads
    fasta_ext = ["fasta", "fa", "fna", "fsa", "fas", "seq"]
    if os.path.isdir(all_reads):
        # this takes effect if the args.reads_dir==True
        rct = 0
        for file in os.scandir(all_reads):
            if file.path.split(".")[-1] in fasta_ext:
                try:
                    rcounter = get_reads_by_header(header, file, non_int_file)
                    rct += sum(rcounter.values())
                except TypeError:
                    ilog.vlprint("Empty read file or unsupported file type", 3)
        ilog.vlprint(f"Retrieved {rct} non-int-reads", 2)
    else:
        rct = get_reads_by_header(header, all_reads, non_int_file)
        ilog.vlprint(f"Retrieved {sum(rct.values())} non-int-reads", 2)

    mm2_out = "confirm_non_int_reads.paf"
    ilog.vlprint("coverage check for evaluation: mapping reads", 5)
    try:
        minimap2_sp(args.genome, non_int_file, mm2_out, args, secondary="no")
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, f"Minimap2 for confirmation of non-int_reads failed",
            exit=False
            )

    df, cols = make_mm_paf_df(mm2_out, ava=False, quick_mm2=False)
    fildf = filt_sort_paf_df(df, cols, mm2_mapq=False, no_culling=False)
    confirmed = defaultdict(list)
    not_conf = defaultdict(list)

    piso_by_name = {p.name:p for p in prim_int_site_objs}
    for mapping in fildf:
        qname = mapping['qname']
        if qname in read2loc:
            chrom = mapping["tname"]
            if chrom in read2loc[qname]:
                if mapping["tstart"]<read2loc[qname][chrom]<mapping["tend"]:
                    for site in read2site[qname]:
                        confirmed[site].append(qname)
            else:
                for site in read2site[qname]:
                    not_conf[site].append(qname)
        else:
            ilog.vlprint(
                f"WARNING:{qname} is not in read2loc confirm_non_int_reads", 5
                )
    
    for iso in piso_by_name:
        # if iso in confirmed:
        piso = piso_by_name[iso]
        piso.conf_non_int_als = confirmed[iso]
        conf = [r for r in piso.pot_non_int_als if r in confirmed[iso]]
        new = [r for r in confirmed[iso] if r not in piso.pot_non_int_als]
        not_conf = [r for r in piso.pot_non_int_als if r not in confirmed[iso]]
        piso.calc_total_cov(mode="confirmed")
        if any([len(new)>0, len(not_conf)>0]):
            ilog.vlprint(
                f"{piso.name} conf:{len(conf)} new:{len(new)} not_conf:"
                f"{len(not_conf)} total_cov: {piso.total_cov}", 6
            )

    # include int_als supporting prox_hom int-site into conf_non_int_als of site 
    # on other homolog
    all_ph_intreads, prox_hom = get_prox_hom_int_alnms(prim_int_site_objs)
    for piso in prox_hom:
        matches = [r for r in all_ph_intreads if r.CHR==piso.chrom]
        matches = [r for r in matches if r.READ not in piso.int_als]
        matches = [
            r for r in matches if r.chr_start_end[0]<piso.approx_loc<r.chr_start_end[1]
            ]
        matches = {r.READ: r for r in matches}
        piso.proxhom_non_int_als.update(matches)
        piso.calc_total_cov(mode="all")
        ilog.vlprint(
                    f"Included {len(matches)} reads from proximal site on homolog into {piso.name} total_cov: {piso.total_cov}", 6
                )


def mmap_ava_comp_primary(outdir, site_objs, primary_int_site_objs, args):
    """Evaluate primary int-sites for similarities or overlap, by mapping   
    spannned sequence contexts against themselves with minimap2."""
    ilog.vlprint(
        'Using minimap2 for all-vs-all comparison of primary int-sites', 3
        )
    os.chdir(outdir)
    span_seqs = os.path.join(outdir, "eval_int_sites_spanseq.que")
    mm_out = os.path.join(outdir, "eval_ava_comp_prim_sites_mm2.paf")
    fp_prim_pis_objs = {}

    try:
        # subprocess.run(
        #         [
        #         args.minimap2,
        #         "-x", f"ava-{args.seq_tec}",
        #         span_seqs,
        #         span_seqs,
        #         "-o", mm_out
        #         ],
        #     check=True, capture_output=True
        #     )
        minimap2_sp(span_seqs, span_seqs, mm_out, args, map_mode="ava")
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "Minimap2 for ava-overlap of primary integration-sites failed",
            exit=False
            )
    
    try:
        mm_df, cols = make_mm_paf_df(mm_out, ava=True)
        filt_mq = mm_df[mm_df["qname"] != mm_df["tname"]]
        filt_mq = filt_mq[filt_mq["map_len"] > 1000]
        filt_mq = filt_mq[filt_mq["matches"] > 1000]
        filt_mq = filt_mq[filt_mq["matches"]/filt_mq["map_len"] > 0.5]
        filt_mq = filt_mq[filt_mq["matches"]/filt_mq["qlen"] > 0.2]
        with open("eval_ava_comp_prim_sites_mm2_filt.paf", "w") as filt:
            filt_csv = filt_mq.to_csv(index=False, sep="\t")
            filt.write(filt_csv)

            filt_csv = filt_csv.splitlines()
            csv_lists = [line.strip().split("\t") for line in filt_csv]
            csv_lists = [line for line in csv_lists if line[0] != "qname"]

            # convert "IntSite_Pos_spanned_context" name to "IntSite_Pos" 
            # nomenclature for use as key in mm2_false_pos to enable following 
            # sorting and filtering steps

            # generate tuples of the two related sites
            mm2_fptp_pairs = {
                ("_".join(ls[0].split("_")[:-2]), "_".join(ls[5].split("_")[:-2])):ls for ls in csv_lists
            }
            # then check which one is likely to be the major and which the stray
            mm2_false_pos0 = {
                min(k, key=lambda x: len(site_objs[x].alignments)):k for k in mm2_fptp_pairs.keys()
            }
            # assign the associated true_pos as value to the false_pos as key
            mm2_false_pos = {k:[s for s in v if s!=k][0] for k,v in mm2_false_pos0.items()}
            # attach the supporting alignments to the major site, but without
            # recalculating location etc (i.e. not performing full merger), for
            # it is likely that those reads are of lesser quality and would
            # negatively impact on location precision.
            primary_int_site_objs = {pis.name:pis for pis in primary_int_site_objs}
            validated_tp = []
            for pot_fp in mm2_false_pos:
                merge = True
                stray = primary_int_site_objs[pot_fp]
                major = primary_int_site_objs[mm2_false_pos[pot_fp]]
                if stray.chrom==major.chrom:
                    if abs(stray.approx_loc-major.approx_loc) < 500000:
                        # overlap of spanned regions
                        mult_evid = []
                        for al in stray.supp_als:
                            mult_evid.append(al.cand_read_obj.multiple_sep_int)
                            mult_evid.append(al.cand_read_obj.multimer_int)
                        for al in major.supp_als:
                            mult_evid.append(al.cand_read_obj.multiple_sep_int)
                            mult_evid.append(al.cand_read_obj.multimer_int)
                            # print("major_evidence_tmp", major_evidence_tmp)
                        if any(mult_evid):
                            # evidence of multiple integrations 
                            # in at least one of both
                            merge=False
                            major.multi_site, stray.multi_site = True, True
                        else:
                            # the threshold value below could either be set by a
                            #  parameter, or be statistically derived by analyzi
                            # ng the deviation of alignment edges used for calcu
                            # lation of integration-site location
                            if len(stray.supp_als)/len(major.supp_als) < 0.5:
                                # the big difference in coverage of the two sites
                                # in absence of evidence for a multimeric/clustered
                                # site indicates a false-positive site induced by 
                                # alignment artifacts/errors
                                merge = True  
                            else:
                                # the two sites likely lie on different 
                                # chromosome homologs
                                merge = False
                                major.prox_homolog.append(stray.name)
                                stray.prox_homolog.append(major.name)
                if merge:
                    major.supp_als.extend(stray.supp_als)
                    major.supp_absorbed_fp.extend(stray.supp_als)
                else:
                    validated_tp.append(pot_fp)

        mm2_false_pos = {k:v for k,v in mm2_false_pos.items() if k not in validated_tp}
        fp_prim_pis_objs = {k:v for k,v in site_objs.items() if k in mm2_false_pos}
        site_objs = {k:v for k,v in site_objs.items() if k not in mm2_false_pos}
        primary_int_site_objs = [
            pis for pis in primary_int_site_objs.values() if pis.name not in mm2_false_pos
            ]
        for fps in mm2_false_pos:
            fp_prim_pis_objs[fps].mm2_false_pos_evidence = {fps: mm2_fptp_pairs[mm2_false_pos0[fps]]}
              
    except pd.errors.EmptyDataError as e:
        ilog.vlprint("minimap_ava_comp_primary did not yield any hits", 2)

    return primary_int_site_objs, fp_prim_pis_objs
        

def mmap_eval_primary(ref_genome, outdir, sites_by_chrom, args):
    """Evaluate primary int-sites by mapping cand_reads (without integration) vs
     ref_genome with minimap2."""
    ilog.vlprint('Using minimap2 to evaluate primary int-sites', 3)
    os.chdir(outdir)
    cand_reads_wo = os.path.join(outdir, "cand_reads_without_integration.fasta")
    mm_out = os.path.join(outdir, "eval_prim_sites_mm2.paf")

    try:
        # subprocess.run(
        #         [
        #         args.minimap2,
        #         "-x", f"map-{args.seq_tec}",
        #         ref_genome,
        #         cand_reads_wo,
        #         "-o", mm_out
        #         ],
        #         check=True, capture_output=True
        #     )
        minimap2_sp(ref_genome, cand_reads_wo, mm_out, args)
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "Minimap2 for evaluation of primary integration-sites failed",
            exit=False
            )
    try:
        mm_df, cols = make_mm_paf_df(mm_out)
        mm_alnm = filt_sort_paf_df(mm_df, cols)

        intinread_pos = {}
        for chrom in sites_by_chrom:
            for site in sites_by_chrom[chrom]:
                for sal in sites_by_chrom[chrom][site].supp_als:
                    intinread_pos[sal.cand_read_obj.RNAME] = sal.cand_read_obj.read_cov_range[0]

        aln_ct = 0
        sort_count = defaultdict(list)
        excluded = [*args.filt_for_len,*args.nosup_passint,*args.no_genseq]
        for aln in mm_alnm:
            aln_ct += 1
            if aln["qname"] not in excluded:
                chrom = aln["tname"]
                int_pos = None
                for loc in sites_by_chrom[chrom]:
                    try:
                        int_pos = intinread_pos[aln["qname"]]
                    except KeyError:
                        ilog.vlprint(
                            f"INFO:{aln['qname']} not in intinread_pos.",6
                            )
                    if int_pos:
                        if aln["strand"]=="+":
                            distloc = abs(aln["tstart"]-aln["qstart"]+int_pos-loc)  
                        else:
                            distloc = abs(aln["tend"]+aln["qstart"]-int_pos-loc)
                        if distloc <= (int_pos/20)+args.min_dist_ints:
                            sites_by_chrom[chrom][loc].minimap_supp[aln["qname"]] = aln
                            sort_count[aln["qname"]].append(sites_by_chrom[chrom][loc].name)
                    else:
                        mindistloc = min(
                            [abs(aln["tstart"] - loc), abs(aln["tend"] - loc)]
                            )
                        if mindistloc <= int(aln["qlen"]):
                            sites_by_chrom[chrom][loc].minimap_supp[aln["qname"]] = aln
                            sort_count[aln["qname"]].append(sites_by_chrom[chrom][loc].name)
        
        sort_al = len(sort_count)
        ilog.vlprint(
            f"{sort_al}/{aln_ct} mm2-alignments were sorted to PotIntSite-objs",
            4
            )
        for al in sort_count:
            if len(sort_count[al]) > 1:
                ilog.vlprint(
                    f"{al} has been sorted to {len(sort_count[al])} different"
                    f" potential Integration sites: {sort_count[al]}.", 2
                )
    except pd.errors.EmptyDataError as e:
        ilog.vlprint("minimap_ava_comp_primary did not yield any hits", 2)

    return sites_by_chrom


def mmap_for_splitreads(outdir, prim_site_objs, intlen, args):
    """Evaluate primary int-sites by mapping cand_reads vs spanning regions with 
    minimap2."""
    ilog.vlprint(
        'Mapping cand_reads (with integration) using minimap2 to generate split-read info'
        , 3
        )
    os.chdir(outdir)
    cand_reads = os.path.join(outdir, "Integration_candidate_reads.fasta")
    span_seqs = os.path.join(outdir, "eval_int_sites_spanseq.que")
    mm_out = os.path.join(outdir, "eval_prim_sites_mm2_splitreads.paf")

    longjoin_bw = int(intlen * 0.80) # 0.9
    chaining_bw = 500 if longjoin_bw > 899 else intlen * 0.5
    try:
        # subprocess.run(
        #         [
        #         args.minimap2,
        #         # restricting long-join bandwidth to length close to int-size
        #         # with -r allows split-read info generation with minimap2, 
        #         # which otherwise often is kept as large insertions (see cigar)
        #         "-r", f"{chaining_bw},{longjoin_bw}",
        #         "-x", f"map-{args.seq_tec}",
        #         span_seqs,
        #         cand_reads,
        #         "-o", mm_out
        #         ],
        #         check=True, capture_output=True
        #     )
        minimap2_sp(
            span_seqs,
            cand_reads,
            mm_out,
            args,
            optargs=["-r", f"{chaining_bw},{longjoin_bw}",]
            )
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "Minimap2 for split-read mapping to primary int-sites failed",
            exit=False
            )
    try:
        mm_df, cols = make_mm_paf_df(mm_out)
        mm_alnm = filt_sort_paf_df(mm_df, cols, no_culling=True)

        sites_by_name = {piso.name:piso for piso in prim_site_objs}
        tname_qname_dct = {
            "_".join(al["tname"].split("_")[:-2]):defaultdict(list) for al in mm_alnm
            }

        for al in mm_alnm:
            tname = "_".join(al["tname"].split("_")[:-2])
            tname_qname_dct[tname][al["qname"]].append(al)
        
        read_ct = 0
        aln_ct = 0
        for tname in tname_qname_dct:
            if tname in sites_by_name:
                sites_by_name[tname].minimap_supp_sr = tname_qname_dct[tname]
                for qname in tname_qname_dct[tname]:
                    read_ct += 1
                    aln_ct += len(tname_qname_dct[tname][qname])
            
        ilog.vlprint(
            f"{aln_ct} alignments of {read_ct} reads were sorted to PotIntSite-objs",
            4
            )
    except pd.errors.EmptyDataError as e:
        ilog.vlprint("minimap_ava_comp_primary did not yield any hits", 2)

    return prim_site_objs


def compare_prim_to_sec(sec_sites: dict, pot_int_site_objs: list):
    """Compare primary to secondary intsite candidates"""
    ilog.vlprint('Comparing primary to secondary intsite candidates', 1)
    cbs_comp :dict= defaultdict(dict)
    pis_dct = {
        pis.name: pis for pis in pot_int_site_objs if not pis.ssbr == "n.d."
        }
    for site in pis_dct:
        piso = pis_dct[site]
        cum_bs = sum([sum(rtr_al.bitscore) for rtr_al in piso.supp_als])
        cbs_comp[site]["primary"] = cum_bs
    for site in sec_sites:
        if site in pis_dct:
            cbs_comp[site]["secondary"] = sec_sites[site].analyze_sec_site()
    for site in pis_dct:
        ilog.vlprint(f"Alignment score ratios {site}", 3, alt_log="site_specs.log")
        # Two alignment score ratios:
        # 1. Similar (secondary) site bitscore ratio (SSBR)
        # 2. Identical matches ratio (IMR)
        if len(cbs_comp[site]) == 2:
            SSBR = cbs_comp[site]["secondary"] / cbs_comp[site]["primary"]
            pis_dct[site].ssbr = SSBR
        elif len(cbs_comp[site]) == 1:
            pis_dct[site].ssbr = "none"
        else:
            pis_dct[site].ssbr = "cbs_comp error"
        IMR = pis_dct[site].ident_match_ratio
        ilog.vlprint(f"IMR {IMR}", 3, alt_log="site_specs.log")
        ilog.vlprint(f"SSBR {pis_dct[site].ssbr}", 3, alt_log="site_specs.log")
        ilog.vlprint(f"blast_supp {len(pis_dct[site].supp_als)}", 3, alt_log="site_specs.log")
        ilog.vlprint(f"minimap_supp {len(pis_dct[site].minimap_supp)}", 3, alt_log="site_specs.log")
    ilog.vlprint('Finished comparing primary to secondary intsite candidates', 1)
    return pis_dct


def mult_ind_int_read_info(
        prim_isite_objs: list,
        site_objs: dict,
        cand_read_objs: list,
        rtr_alnms: list,
        args: object
        ):
    """Use information from reads containing evidence of multiple independent integrations"""
    multi_cro = {cro.RNAME:cro for cro in cand_read_objs if cro.multiple_sep_int}
    multi_rtr = {rtr.READ:rtr for rtr in rtr_alnms if rtr.READ in multi_cro}
    # single_cro = {cro.RNAME:cro for cro in cand_read_objs if not cro.multiple_sep_int}
    # single_rtr = {rtr.READ:rtr for rtr in rtr_alnms if rtr.READ in single_cro}

    piso_by_chrom :dict= defaultdict(dict)
    for site in prim_isite_objs:
        piso_by_chrom[site.chrom][site.approx_loc] = site

    if len(multi_cro) > 0:
        il_map_figs.gen_multi_int_read_figs(multi_cro)
        read_site_assoc = defaultdict(list)
        for pio in prim_isite_objs:
            for al in pio.supp_als:
                read_site_assoc[al.READ].append(pio.name)
    # site_linkage = defaultdict(set)
    for chrom in piso_by_chrom:
        chrom_multi_rtr = {k:v for k,v in multi_rtr.items() if v.CHR==chrom}
        # site_link_chrom = {}
        allsites_chr = list(piso_by_chrom[chrom].keys())
        for rtr in chrom_multi_rtr:
            alrdy_ass = [int(p.split("_")[-1]) for p in read_site_assoc[rtr]]
            int_pos = sorted(
                multi_cro[rtr].orig_alignment,
                key=lambda x: x[0]
            )
            int_ct = len(int_pos)
            ir = list(range(int_ct-1))
            int_dists = [int_pos[i+1][0]-int_pos[i][1] for i in ir]
            rtr_chrom_pos = sorted(
                chrom_multi_rtr[rtr].read_start_end_chr_start_end,
                key=lambda x: x[0][0]
                )
            # print("chrom", chrom)
            # print(rtr, "read_start_end_chr_start_end", rtr_chrom_pos)
            # print(rtr, "int_pos", int_pos)
            # print("int_dists in read", int_dists, chrom_multi_rtr[rtr].sstrand)
            # print(rtr, chrom_multi_rtr[rtr].ins_gap)
            # print("##############################")
            segments = copy(int_pos)
            rtr_ct = 0
            ins_ct = 0
            for rtr_range in rtr_chrom_pos:
                rtr_ct += 1
                read_range = rtr_range[0]
                ip_ct = 0
                for ip in int_pos:
                    if read_range[1]-ip[0] <= 3:
                        if read_range[0] < ip[0]:
                            break
                        else:
                            ip_ct += 1
                    else:
                        ip_ct += 1
                
                if ip_ct==0:
                    segments.insert(ip_ct, rtr_range)
                    ins_ct += 1
                else:
                    segments.insert(ip_ct+ins_ct, rtr_range)
                    ins_ct += 1
            chrom_multi_rtr[rtr].multi_info["segments"] = segments
            # print("segments", segments)

            idx = -1
            ip_ct = 0
            ins_chrom_pos = {}
            for seg in segments:
                idx += 1
                if type(seg)==list:
                    # int_pos segment
                    ip_ct += 1
                    refstrand = chrom_multi_rtr[rtr].sstrand[0]
                    if refstrand=="minus":
                        idx_left, idx_right = 0, 1
                    else:
                        idx_left, idx_right = 1, 0
                    if idx==0:
                        if type(segments[1])==tuple:
                            ins_chrom_pos[idx] = segments[1][1][idx_right]
                    elif idx==len(segments)-1:
                        if type(segments[-2])==tuple:
                            ins_chrom_pos[idx] = segments[-2][1][idx_left]
                    else:
                        adj_borders = []
                        if type(segments[idx-1])==tuple:
                            adj_borders.append(segments[idx-1][1][idx_left])
                        if type(segments[idx+1])==tuple:
                            adj_borders.append(segments[idx+1][1][idx_right])
                        ins_chrom_pos[idx] = sum(adj_borders)/len(adj_borders)          
                else:
                    # rtr_range segment
                    pass
            chrom_multi_rtr[rtr].multi_info["ins_chrom_pos"] = ins_chrom_pos

            delta_sites = {}
            for ins in ins_chrom_pos:
                ip = ins_chrom_pos[ins]
                delta_tmp = {(ins,s):abs(s-ip) for s in allsites_chr}
                best = min(delta_tmp, key=lambda x: delta_tmp[x])
                delta_sites.update({best:delta_tmp[best]})
            for site in delta_sites:
                idx, site = site
                piso = piso_by_chrom[chrom][site]
                # filter linked int-sites based on the accuracy of accordance
                linked_pos = [k[1] for k,v in delta_sites.items() if v<=args.min_dist_ints]
                # remove self from the list of linked positions
                linked_pos = [k[1] for k in delta_sites.keys() if k[1]!=site]
                # get the respective PotIntSite objects
                linked_piso = [piso_by_chrom[chrom][pos] for pos in linked_pos]
                if site in alrdy_ass:
                    piso.multread_linked.update({p.name:p for p in linked_piso})
                else:
                    piso.mrli_supp[rtr] = [chrom_multi_rtr[rtr], idx]
                    piso.supp_als.append(chrom_multi_rtr[rtr])
                    piso.multread_linked.update({p.name:p for p in linked_piso})
            chrom_multi_rtr[rtr].multi_info["delta_sites"] = delta_sites
    for piso in prim_isite_objs:
        old_name = piso.name
        piso.re_calc_attributes(rename=True)
        new_name = piso.name
        piso.multread_linked = {v.name:v for v in piso.multread_linked.values()}
        if old_name != new_name:
            ilog.vlprint(
                        f"{old_name} was renamed to {new_name} based on recalcu"
                        f"lation of int-position after inclusion of multi-int-"
                        f"read information", 3
                        )
            piso.old_name = old_name
            piso.dump_reads()
            os.remove(f"{old_name}_supporting_reads.fasta")
            site_objs[old_name].name, site_objs[old_name].old_name = new_name, old_name
            site_objs[old_name].assoc_primary_site = new_name
    # to take effect also in piso_by_chrom and site_objs dcts, keys must be reset
    # at the same time mark as multi-site if applicable
    # multi_rep = open("multi_site_report.csv", "w")
    # proxtrans_rep = open("proximal_trans_sites_report.csv", "w")
    # multi_rep.write(f"site,linked cis-sites,linking-reads,multimer-reads")
    site_objs = {v.name:v for v in site_objs.values()}
    for chrom in piso_by_chrom:
        piso_by_chrom[chrom] = {
            v.approx_loc:v for v in piso_by_chrom[chrom].values()
            }

    return site_objs, piso_by_chrom


def write_multi_reports(piso_by_chrom):
    # mark as multi-site if applicable
    multi_rep = open("multi_site_report.csv", "w")
    multi_rep.write(f"site,linked cis-sites,linking-reads,multimer-reads\n")
    proxtrans_rep = open("proximal_trans_sites_report.csv", "w")
    proxtrans_rep.write(f"site,confirmed trans-sites\n")

    for chrom in piso_by_chrom:
        for key in piso_by_chrom[chrom]:
            pis = piso_by_chrom[chrom][key]
            multi_sep_supp = set()
            multimer_supp = set()
            for rtr in pis.mrli_supp:
                # structure of mrli_supp entries: [chrom_multi_rtr[rtr], idx]
                multi_sep_supp.add(pis.mrli_supp[rtr][0])
            for rtr in pis.supp_als:
                if rtr.cand_read_obj.multiple_sep_int:
                    multi_sep_supp.add(rtr)
                if rtr.cand_read_obj.multimer_int:
                    multimer_supp.add(rtr)
            pis.multi_sep_supp = multi_sep_supp
            pis.multimer_supp = multimer_supp
            if any([any(multi_sep_supp), any(multimer_supp)]):
                pis.multi_site = True
                # names of recalculated linked sites occuring later in the loop
                # in mult_ind_int_read_info, will be outdated; thus reset here
                pis.multread_linked = {v.name:v for v in pis.multread_linked.values()}
            if pis.multi_site:
                multi_rep.write(
                    f"{pis.name},{';'.join(pis.multread_linked.keys())},"
                    f"{len(pis.multi_sep_supp)},{len(pis.multimer_supp)}\n"
                    )
            if len(pis.prox_homolog) > 0:
                proxtrans_rep.write(f"{pis.name},{';'.join(pis.prox_homolog)}\n")
    multi_rep.close()
    proxtrans_rep.close()


def write_loss_gain_reports(prim_isite_objs):
    entries = 0
    with open("seq_gain_loss_at_intsites.csv", "w") as rep:
        rep.write(
            "Integration,seq loss [min/median/most_freq](bp), seq gain left [min/median"\
            "/most_freq] (bp), seq gain right [min/median/most_freq] (bp), confounding factors\n"
        )
        for p in prim_isite_objs:
            p.calc_del_at_junction()
            p.calc_seq_gain()
            mls = p.seq_loss_info["min_loss"] 
            medls = p.seq_loss_info["med_loss"]
            mfls = p.seq_loss_info["most_freq_gap"]
            gl = p.seq_gain_info["left"]
            mgl, medgl, mfgl = gl["min"], gl["median"], gl["most_freq"]
            gr = p.seq_gain_info["right"]
            mgr, medgr, mfgr = gr["min"], gr["median"], gr["most_freq"]
            all_vals = [mls, medls, mfls, mgl, medgl, mfgl, mgr, medgr, mfgr]
            if mfls==0 and medls==1:
                all_vals[1] = 0
            if mfgl==0 and medgl==1:
                all_vals[4] = 0
            if mfgr==0 and medgr==1:
                all_vals[-2] = 0
            all_vals = [v for v in all_vals if type(v)!=str]
            pos_vals = [v for v in all_vals if v>=1]
            confound = "multi_site" if p.multi_site else "None"
            if any([pos_vals]):
                entries += 1
                rep.write(
                    f"{p.name},{mls}/{medls}/{mfls},{mgl}/{medgl}/{mfgl},"
                    f"{mgr}/{medgr}/{mfgr},{confound}\n"
                    )
    if not entries:
        os.remove("seq_gain_loss_at_intsites.csv")
    

def run_il_evaluate(
    outdir: str,
    ref_genome: str,
    prim_isite_objs: list,
    cand_read_objs: list,
    rtr_alnms: list,
    all_reads: str,
    intlen: int,
    thr: int,
    args
    ):
    """Execution of evaluation functions"""
    ilog.verbosity = args.verbosity
    ilog.vlprint("Starting integration-site evaluation module", 5)
    report_time("evaluation module")
    alignments = load_alignments(outdir)
    site_objs = read_out_spanned_regions(alignments, prim_isite_objs)
    if not args.no_multi_ints:
        ilog.vlprint("Analyzing multi-int reads", 1)
        site_objs, piso_by_chrom = mult_ind_int_read_info(
            prim_isite_objs, site_objs, cand_read_objs, rtr_alnms, args
            )
    else:
        piso_by_chrom = defaultdict(dict)
        for site in prim_isite_objs:
            piso_by_chrom[site.chrom][site.approx_loc] = site
        
    mmap_eval_primary(ref_genome, outdir, piso_by_chrom, args)
    spanned_seqs = get_spanned_seqs(site_objs, ref_genome)
    simil_search_res = search_similar_refgen(spanned_seqs, ref_genome, args)
    prim_isite_objs, mm2_fp_pis_objs = mmap_ava_comp_primary(
                                                        outdir,
                                                        site_objs,
                                                        prim_isite_objs,
                                                        args
                                                        )
    ilog.vlprint(
    f"{len(mm2_fp_pis_objs)} potential integration site candidate/s"
    f" was/were determined to be false positives due to extensive overlap"
    f" and similarity.", 1
    )
    site_objs = filt_similar_seqs(site_objs, simil_search_res)
    check_similar_regions(site_objs)
    sec_site_objs = create_alt_reg_objs(site_objs, ref_genome)
    blast_cand_reads_vs_sec_site_reg(sec_site_objs, args)
    chk_cr_vs_ssr_success(prim_isite_objs)
    compare_prim_to_sec(sec_site_objs, prim_isite_objs)
    ilog.vlprint(
        f"{len(prim_isite_objs)} integration sites have been confirmed "
        f"after evaluation.", 1
        )
    try:
        mmap_total_cov_primsite(all_reads, prim_isite_objs, args) # site_objs, 
    except Exception as e:
        ilog.vlprint(e, 2)
        ilog.vlprint(
            f"WARNING: Determination of total coverage at in-sites failed.", 2
            )
    confirm_non_int_reads(prim_isite_objs, all_reads, args)
    try:
        mmap_for_splitreads(outdir, prim_isite_objs, intlen, args)
        il_map_figs.gen_map_figs(prim_isite_objs, aligner="minimap2")
    except Exception as e:
        ilog.vlprint(e, 2)
        ilog.vlprint(
            f"WARNING: Creation of figs for candidate-read mappings failed.", 2
            )
    if not args.no_multi_ints:
        write_multi_reports(piso_by_chrom)
    try:
        write_loss_gain_reports(prim_isite_objs)
    except:
        ilog.vlprint(
            "Report for seq. loss/gain at int. sites could not be generated", 1
            )

    report_time("Evaluation module", end=True)
    return prim_isite_objs
    

if __name__ == "__main__":
    outdir = sys.argv[1]
    ref_genome = sys.argv[2]
    pickled_primary_intsites = sys.argv[3]
    # run_il_evaluate(outdir, ref_genome, pickled_primary_intsites)
  