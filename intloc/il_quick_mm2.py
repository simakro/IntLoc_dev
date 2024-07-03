# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

from collections import defaultdict
import subprocess
import os
import pandas as pd

from il_logging import Ilogger
from il_chk_conv import report_time, get_reads_by_header
from il_report import write_int_csv, draw_generic_ideogram
from il_evaluate import make_mm_paf_df, filt_sort_paf_df
from il_subprocesses import caplog_sp_error, minimap2_sp
# from integration_locator import BlastCandidateRead, PotIntSite
from il_classes import MmapCandidateRead, PotIntSite
from integration_locator import write_intread_summary_txt, report_reads_info
from integration_locator import merge_pot_int_sites, contains_multimers
from integration_locator import write_reads_wo_integration, filter_potintsites


ilog = Ilogger()
ilog.module = __name__

def mm2_cand_reads(construct, reads, outdir, args, read_dir=False, bait=False):
    """Identify candidate reads containing integration of interest in sequencing
     data"""
    ilog.vlprint(
        "Identifying candidate reads containing integration of interest in"
        " sequencing data using minimap2", 0
        )
    readpath = os.path.split(reads)
    workpath = readpath[0]
    os.chmod(workpath, 0o777)
    os.chdir(outdir)

    if all([args.paftools, args.samtools]):
        presets = "-ax"
        filetype = "sam"
    else:
        presets = "-x"
        filetype = "paf"

    if bait:
        # there could be downstream compatibility problems with bait.paf files; check!
        cand_reads = os.path.join(outdir, f"bait_candidate_reads.{filetype}")
    else:
        cand_reads = os.path.join(outdir, f"candidate_reads.{filetype}")

    if os.path.exists(cand_reads) and not read_dir:
        ilog.vlprint(
            'candidate_reads file exists. '
            'Skipping search for reads containing construct', 1
            )
    else:    
        ilog.vlprint(
            'Starting minimap2-search to identify reads containing '
            'integration-construct sequences', 1
            )
        try:
            # subprocess.run(
            #     [
            #     args.minimap2,
            #     presets, f"map-{args.seq_tec}",
            #     "--secondary=no", # available in minimap2 2.3 and higher
            #     "--sam-hit-only",
            #     construct,
            #     reads,
            #     "-o", cand_reads  # available in minimap2 2.16 and higher
            #     ],
            #     check=True, capture_output=True
            # )
            minimap2_sp(
                construct,
                reads,
                cand_reads,
                args,
                secondary="no",
                preset=presets,
                optargs=["--sam-hit-only"]
                )
        except subprocess.CalledProcessError as e:
            caplog_sp_error(
                e, "minimap2-search for identification of candidate reads failed"
                )
        except FileNotFoundError as fe:
            ilog.vlprint(f"mm2_cand_reads: FileNotFoundError: {fe}", 1)
    with open(cand_reads, "r") as crf:
        data = crf.read().split("\n")
    ilog.vlprint(f"Obtained {len(data)} hits in candidate read search", 0)
    return cand_reads


def extract_candidate_reads(cand_reads):
    """Extract candidate reads from sam-file using samtools"""
    ilog.vlprint(
                'Extracting Integration candidate reads using samtools', 0
                )
    try:
        sbp_out = subprocess.run(
            [
            "samtools",
            "fasta",
            cand_reads,
            "-F", "4,256,512,2048"
            ],
            check=True, capture_output=True
        )
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "extraction of candidate reads with samtools failed"
            )
    except FileNotFoundError as fe:
        ilog.vlprint(f"extract_candidate_reads: FileNotFoundError: {fe}", 1)

    with open("Integration_candidate_reads.fasta", "w") as out:
        cr_data = sbp_out.stdout.decode()
        out.write(cr_data)
        found_ct = [n for n in range(int(len(cr_data.split("\n"))/2))]
    return found_ct


def get_qnames_from_sam_or_paf(cand_reads):
    """Extracts query-names from sam or paf alignment file"""
    ilog.vlprint(
        f"Get query-names from sam or paf alignment file for read extraction"
        f" using custom intloc functions", 2
    )
    candread_names = set()
    with open(cand_reads, "r") as sam:
        for line in sam:
            if line.startswith(("@", "#","qname")):
                pass
            else:
                ls = line.strip().split("\t")
                if len(ls) > 0:
                    candread_names.add(ls[0])
    return list(candread_names)
        

def convert_sam2paf(sam_file):
    stem = "".join(sam_file.split(".")[:-1])
    ilog.vlprint(
        f'Convert {os.path.split(stem)[-1]} reads sam to paf using paftools', 0
        )
    success = False
    try:
        sbp_out = subprocess.run(
            [
            "paftools.js",
            "sam2paf",
            sam_file,
            ],
            check=True, capture_output=True
        )
        success = True
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, f"Conversion of {stem} sam to paf with paftools failed"
            )
    except FileNotFoundError as fe:
        ilog.vlprint(f"convert_sam2paf: FileNotFoundError: {fe}", 1)

    out_file = f"{stem}.paf"
    with open(out_file, "w") as out:
        out.write(sbp_out.stdout.decode())

    return  out_file


# class MmapCandidateRead(BlastCandidateRead):
#     """Class for creation of minimap2 candidate read objects"""
#     def __init__(
#         self, ID, RNAME, slen, length, constr_cov_range, read_cov_range, pident,  
#         mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus, constr_len, seq,
#         orig_alignment, strand, matches, mapq, al_type, minimizer, chain_score_prim,
#         seq_div, len_rep_seeds
#         ):
#         super().__init__(
#             ID, RNAME, slen, length, constr_cov_range, read_cov_range,   
#             pident, mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus, 
#             constr_len, seq, orig_alignment, aligner="minimap2"
#         )
        
#         self.strand = strand
#         self.mapq = mapq
#         self.al_type = al_type
#         self.minimizer = minimizer
#         self.chain_score_prim = chain_score_prim
#         self.seq_div = seq_div
#         self.len_rep_seeds = len_rep_seeds


def gen_mmap_candread_inst(cand_dcts, args):
    cand_dcts = {dct["qname"]: dct for dct in cand_dcts}
    crs_by_name = {}
    with open("Integration_candidate_reads.fasta", "r") as icrs:
        data = icrs.read().split(">")
        data.pop(0)
        for header_seq in data:
            header, seq = header_seq.strip().split("\n")
            crs_by_name[header] = seq.strip()
            try:
                cand_dcts[header]["seq"] = seq
            except:
                continue
    
    mm_cris = []
    id_ct = 0
    for crd in cand_dcts:
            crd = cand_dcts[crd]
            id_ct += 0
            new_inst = MmapCandidateRead(
                                        id_ct,
                                        crd["qname"],
                                        crd["qlen"],
                                        crd["map_len"],
                                        (int(crd["tstart"]), int(crd["tend"])),
                                        (int(crd["qstart"]), int(crd["qend"])),
                                        False, # place for pident
                                        False, False, False, False, False, False, 
                                        # mismatch,gapopen,evalue,bitscore,qcovhsp,qcovus
                                        crd["tlen"],
                                        crd["seq"],
                                        list(), # False, #orig_alignment
                                        crd["strand"],
                                        crd["matches"],
                                        crd["MAPQ"],
                                        crd["al_type"],
                                        crd["minimizer"],
                                        crd["chain_score_prim"],
                                        crd["seq_div"],
                                        crd["len_rep_seeds"],
                                        args
                                    )
            mm_cris.append(new_inst)
    return mm_cris


def mmquick_int_genome(genome, outdir, args):
    """Identify potential integration locations by mapping  
    candidate reads against the reference using minimap2."""

    ilog.vlprint(
        'Starting minimap2-search to identify location of '
        'integrations in reference genome', 1
        )
    os.chdir(outdir)
    if args.no_culling:
        int_cand_reads = "Integration_candidate_reads.fasta"
    else:
        int_cand_reads = "cand_reads_without_integration.fasta"

    outfile = os.path.join(outdir, "Integration_sites.paf")
    try:
        # subprocess.run(
        #         [
        #         args.minimap2,
        #         "--secondary=no", # available in minimap2 2.3 and higher
        #         # "-N", "0", # can reduce mapping accuracy (def=5) use --secondary
        #         "-x", f"map-{args.seq_tec}",
        #         genome,
        #         int_cand_reads,
        #         "-o", outfile # available in minimap2 2.16 and higher
        #         ],
        #     check=True, capture_output=True
        #     )
        minimap2_sp(genome, int_cand_reads, outfile, args, secondary="no")
        with open(outfile, "r") as out:
            data = out.read().split("\n")
            refhit_mappings = [l.split("\t")[0] for l in data if len(l)>0]
            refhit_reads = set(refhit_mappings)
        ilog.vlprint(
            f"Generated {len(refhit_mappings)} reference mappings of"
            f" {len(refhit_reads)} integration reads with minimap2.", 0
            )
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "Minimap2 for int-sites identification failed",
            exit=False
            )
    except FileNotFoundError as fe:
        ilog.vlprint(f"mmquick_int_genome: FileNotFoundError: {fe}", 1)
    
    return outfile


def filt_mm_rtr_als(in_paf, args):
    try:
        mm_df, cols = make_mm_paf_df(in_paf)
        mm_alnm = filt_sort_paf_df(mm_df, cols, mm2_mapq=0, no_culling=args.no_culling)
        rhr_after_filt = set([aln["qname"] for aln in mm_alnm])
        ilog.vlprint(
            f"{len(rhr_after_filt)} mm2-aligned reads remain after filtering",4
            )
        return mm_alnm
    except pd.errors.EmptyDataError as e:
        ilog.vlprint(
            "mmap_int_genome for identification of int-sites has failed", 0
            )
        return []


def amend_mm_read2ref_als(mm_alnm, mm_candread_ins):
    """Attach cand_read_objects to rtr-alignments and calculate approx.
     insertion positions in chromosom"""
    ilog.vlprint(
        "Attaching cand_read_objects to rtr-alignments and calculate approx."
        " insertion positions in chromosom", 0
        )
    mm_alnm = {al["qname"]:al for al in mm_alnm}
    not_in_r2r_als = []
    for cro in mm_candread_ins:
        try:
            pos_in_read_wo = cro.read_cov_range[0]
            # chrom_name = mm_alnm[cro.RNAME]["tname"]
            tstart = mm_alnm[cro.RNAME]["tstart"]
            tend = mm_alnm[cro.RNAME]["tend"]
            qstart = mm_alnm[cro.RNAME]["qstart"]
            qend = mm_alnm[cro.RNAME]["qend"]
            qlen = mm_alnm[cro.RNAME]["qlen"]

            # pos_in_chrom = tstart + (pos_in_read_wo - qstart)

            # if mm_alnm[cro.RNAME]["strand"] == "+":
            #     pos_in_chrom = tstart + (pos_in_read_wo - qstart)
            # else:
            #     pos_in_chrom = tend + (pos_in_read_wo - qstart)

            if mm_alnm[cro.RNAME]["strand"] == "+":
                if not cro.trunc_constr:
                    pos_in_chrom = tstart + (pos_in_read_wo - qstart) # 53-215 bp off 
                    # pos_in_chrom = tstart + (pos_in_read_wo) # schlechter 90-290 bp off
                    # pos_in_chrom = tstart + (qlen - pos_in_read_wo - qstart) # much worse
                elif cro.trunc_constr == 5:
                    pos_in_chrom = tstart + (pos_in_read_wo - qstart) # 0-25bp off
                elif cro.trunc_constr == 3:
                    pos_in_chrom = tend + (qlen - qend) # 0-22bp off
                else:
                    pos_in_chrom = tstart + (pos_in_read_wo - qstart) 
            else:
                if not cro.trunc_constr:
                    # pos_in_chrom = tstart + (qlen - pos_in_read_wo) # around 200bp off
                    # pos_in_chrom = tstart + (pos_in_read_wo - qstart) # very bad 3000-4000 bp off
                    # pos_in_chrom = tstart + (qlen - pos_in_read_wo - qstart) # around 160bp off
                    # pos_in_chrom = tend + (qlen - pos_in_read_wo - qstart) # very bad 12000 bp off
                    # pos_in_chrom = tend - (qlen - pos_in_read_wo - qstart) # very bad 3500 bp off
                    pos_in_chrom = tend - (pos_in_read_wo - qstart) # 106-293 bp off
                    # pos_in_chrom = tend - (pos_in_read_wo) # about 130bp off
                    # pos_in_chrom = tend - (pos_in_read_wo + qstart) # about 170bp off
                elif cro.trunc_constr == 5:
                    # pos_in_chrom = tend + (pos_in_read_wo - qstart) # around 50-150 off
                    pos_in_chrom = tend + (pos_in_read_wo) # around 28-112 off; much better than above
                elif cro.trunc_constr == 3:
                    pos_in_chrom = tstart + (pos_in_read_wo - qlen) # 29-69bp off
                    # pos_in_chrom = tstart + (pos_in_read_wo - qend) # schlechter 50-150bp off
                else:
                    pos_in_chrom = tend + (pos_in_read_wo - qstart)

            mm_alnm[cro.RNAME]["cand_read_obj"] = cro
            mm_alnm[cro.RNAME]["int_pos_chrom"] = pos_in_chrom
        except KeyError:
            not_in_r2r_als.append(cro)
            ilog.vlprint(
                f"{cro.RNAME} is not present in mm2-r2r alnms", 4, logging=False
                )
    ilog.vlprint(
        f"{len(not_in_r2r_als)} cand-reads were not present in mm2-r2r alnms", 4
        )
    return mm_alnm.values()


def mmquick_gen_int_sites(rtr_lst, args):
    """Determine Integration sites from mm2 read to reference alignments"""
    ilog.vlprint(
        "Determining Integration sites from mm2 read to reference alignments", 1
    )
    chr_dict = defaultdict(list)
    loc_dict = {}
    for al in rtr_lst:      
            chr_dict[al["tname"]].append(al)

    for chrom in chr_dict:
        loc_dict[chrom] = {}
        for pot_int_site in chr_dict[chrom]:
            inserted = False
            for pos in loc_dict[chrom]:
                new_pis = pot_int_site["int_pos_chrom"]
                distance = abs(new_pis - pos)
                if distance < args.min_dist_ints:
                    loc_dict[chrom][pos].append(pot_int_site)
                    inserted = True
                    break
            if not inserted:
                loc_dict[chrom][pot_int_site["int_pos_chrom"]] = [pot_int_site]
    ilog.vlprint(f'loc_dict: {loc_dict}', 2, logging=False)

    # output approx. insertion locations of constr. for each read
    for chrom in loc_dict:
        positions = {}
        pos_full_minus = []
        pos_full_plus = []
        pos_5_minus = []
        pos_5_plus = []
        pos_3_minus = []
        pos_3_plus = []

        for pos in loc_dict[chrom]:
            for al in loc_dict[chrom][pos]:
                positions[al["qname"]] = al["int_pos_chrom"]
                if al["strand"]=="-":
                    if not al["cand_read_obj"].trunc_constr:
                        pos_full_minus.append(al["int_pos_chrom"])
                    elif al["cand_read_obj"].trunc_constr == 5:
                        pos_5_minus.append(al["int_pos_chrom"])
                    elif al["cand_read_obj"].trunc_constr == 3:
                        pos_3_minus.append(al["int_pos_chrom"])
                    else:
                        pass
                else:
                    if not al["cand_read_obj"].trunc_constr:
                        pos_full_plus.append(al["int_pos_chrom"])
                    elif al["cand_read_obj"].trunc_constr == 5:
                        pos_5_plus.append(al["int_pos_chrom"])
                    elif al["cand_read_obj"].trunc_constr == 3:
                        pos_3_plus.append(al["int_pos_chrom"])
                    else:
                        pass
   
        positions = sorted(positions.items(), key=lambda x: x[1])
        pos = [x[1] for x in positions]
        reads = [x[0] for x in positions]

    # test for read instances in loc_dict
    pot_int_site_objs = list()
    for chrom in loc_dict:
        for pot_int_site in loc_dict[chrom]:
            # use only the ReadToRefAl instances (tup[1]) instead of the whole  
            # tuple, since the coordinate and insertion-gap inormation is  
            # accesible in these objects via the self.int_coord and 
            # self.ins_gap attributes
            supp_als = {al["qname"]:al for al in loc_dict[chrom][pot_int_site]}
            supp_reads = {al["qname"]:al["cand_read_obj"] for al in supp_als.values()}
            for al in supp_als:
                supp_als[al].pop("cand_read_obj")
            al_coord = [al["int_pos_chrom"] for al in supp_als.values()]
            chrom_coord = (min([al_coord]), max(al_coord))
            new_ins = PotIntSite(
                                chrom,
                                chrom_coord,
                                [],
                                False,
                                args
                                )
            new_ins.minimap_supp.update(supp_als)
            new_ins.mm_only_cand_reads.update(supp_reads)
            new_ins.calc_mm2_loc()
            new_ins.re_calc_attributes()
            new_ins.rename()
            ilog.vlprint(
                f"Generated mm-only PotIntSite {new_ins.name} with"
                f" {len(new_ins.total_support)} supporting reads", 
                3, alt_log="site_specs.log"
                )
            pot_int_site_objs.append(new_ins)
    
    ilog.vlprint(
        f"Generated {len(pot_int_site_objs)} potential integration site objects"
        , 3
    )

    return pot_int_site_objs


def run_quick_mm2(reads, construct, genome, outdir, min_supp, args):
    start = report_time("il_quick_mm2", end=False)
    cand_reads = mm2_cand_reads(construct, reads, outdir, args, read_dir=False, bait=False)

    # try:
    #     found_ct = extract_candidate_reads(cand_reads)
    # except:
    #     search_list = get_qnames_from_sam_or_paf(cand_reads)
    #     found_ct = get_reads_by_header(search_list, reads, "Integration_candidate_reads.fasta")
    # try:
    #     cand_paf = convert_sam2paf(cand_reads)
    # except:
    #     cand_paf = custom_sam2paf_conversion(cand_reads)

    if cand_reads.split(".")[-1]=="sam":
        found_ct = extract_candidate_reads(cand_reads)
        cand_paf = convert_sam2paf(cand_reads)  
    else:
        cand_paf = cand_reads
        search_list = get_qnames_from_sam_or_paf(cand_reads)
        crf_counter = get_reads_by_header(
            search_list, args.reads_file, "Integration_candidate_reads.fasta"
            )
        found_ct = [k for k in crf_counter if crf_counter[k]>0]

    ilog.vlprint(f"Extracted {len(found_ct)} integration candidate reads", 0)
    cand_df, cols = make_mm_paf_df(cand_paf, quick_mm2=True)
    cand_dcts = filt_sort_paf_df(cand_df, cols, mm2_mapq=False, no_culling=False)
    with open("candidate_reads_filt.tmp", "w") as tmp:
        for dct in cand_dcts:
            tmp.write(f"{dct['qname']}\t{dct['qlen']}\n")
    mm_candread_ins = gen_mmap_candread_inst(cand_dcts, args)
    write_reads_wo_integration(mm_candread_ins)
    full_int_len, trunc_constr_len, prov_constr_len = write_intread_summary_txt(mm_candread_ins, args)
    int_len = full_int_len if len(full_int_len) > 0 else trunc_constr_len
    if len(int_len) > 0:
        int_len = sum(int_len)/len(int_len)
    else:
        int_len = 0
    # read2ref_paf = mmap_reads_to_ref("Integration_sites", int_len, args)
    # mmap_int_genome(genome, outdir, pot_int_sites, args)
    read2ref_paf = mmquick_int_genome(genome, outdir, args)
    filt_rtr_als = filt_mm_rtr_als(read2ref_paf, args)
    # rtr_df, cols = make_mm_paf_df(read2ref_paf, quick_mm2=False) #Attention here with cols and quick option something may be off here
    # rtr_dcts = filt_sort_paf_df(rtr_df, cols, mm2_mapq=False, no_culling=False)
    amended_rtr_als = amend_mm_read2ref_als(filt_rtr_als, mm_candread_ins)
    reads_info = {}
    mmer_dct = contains_multimers(mm_candread_ins)[1]
    av_int_len, av_crl = report_reads_info(reads_info, found_ct, mm_candread_ins, mmer_dct, full_int_len, trunc_constr_len, prov_constr_len, args)
    # pot_int_objs = gen_mmonly_integration_sites(amended_rtr_als, mm_candread_ins, args)
    pot_int_objs = mmquick_gen_int_sites(amended_rtr_als, args)
    ilog.vlprint(f"Generated {len(pot_int_objs)} potential integration sites",1)
    perf_merge = True
    iterations = 0
    while perf_merge and iterations < 10:
        iterations += 1
        pot_int_objs, perf_merge = merge_pot_int_sites(pot_int_objs)
    ilog.vlprint(f"{len(pot_int_objs)} pot.-int.-sites remain after merging",1)
    pot_int_objs = filter_potintsites(pot_int_objs, min_supp, reads_info, args)
    ilog.vlprint(f"{len(pot_int_objs)} PotIntSites remain after filtering",1)
    # write_output_txt(pot_int_objs)
    for ints in pot_int_objs:
        # ints.supp_als = list(ints.minimap_supp.values())
        ilog.vlprint(f"{ints.name} mm2_supp: {len(ints.supp_als)}",2)
    int_dict = write_int_csv(pot_int_objs)
    # if not args.species:
    draw_generic_ideogram(int_dict, pot_int_objs, args)
    end = report_time("il_quick_mm2", end=True)
    intcnt = len(pot_int_objs)
    multmers = None
    # Detection of multimers not attempted: the cand-read-method
    # check_for_concatemers() has not been tested with MmapCandidateRead yet,
    # although it has been inhereited from BlastCandidateRead and theoretically 
    # could/should work; currently Evidence for multimers is not reported
    # for il_quick_mm2 and il_intra

    return pot_int_objs, av_crl, start, end, intcnt, multmers, reads_info