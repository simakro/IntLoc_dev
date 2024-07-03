# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
import sys
import json
import subprocess
from collections import defaultdict

from il_logging import Ilogger
from il_dumps import write_reads_wo_integration, generate_spec_ideogram
from il_subprocesses import caplog_sp_error, minimap2_sp
from il_evaluate import make_mm_paf_df, filt_sort_paf_df
from il_chk_conv import get_full_species_name, get_reads_by_header
from il_report import write_int_csv, draw_generic_ideogram
from il_quick_mm2 import *
from il_classes import PotIntSite
from il_aux.il_integrator import Contig
from integration_locator import get_reads_by_header, contains_multimers


ilog = Ilogger()
ilog.module = __name__ 

# 1. minimap construct/intra-seq to ref-genome to find all occurences (save locations in file)
# 2. generate a genome depleted/masked for the construct insertions
# 3. map sequencing reads to the genome and identify all integrations
# 4. compare all found integrations to list of known integrations (generated in step 1)
# 5. the difference contains all new integratrions => report both
##
# The only abundant LINE in humans is LINE1. The human genome contains an estimated 
# 100,000 truncated and 4,000 full-length LINE-1 elements.

def find_ref_intra_sites(genome, args, checking=False):
    """Identify pre-existing integrations of the query sequence in the genome"""
    ilog.vlprint(
        "Identifying pre-existing integrations of the query sequence in the"
        " genome using minimap2", 0
        )
    if not checking:
        out_file = "intra_sites.paf"
    else:
        out_file = "check_masking.paf"

    constr_path = os.path.split(args.construct)[0]
    os.chmod(constr_path, 0o777)
    os.chdir(args.outdir)
    intra_sites = os.path.join(args.outdir, out_file)

    try:
        # subprocess.run(
        #     [
        #     args.minimap2,
        #     "-x", f"map-{args.seq_tec}",
        #     "--secondary=no", # available in minimap2 2.3 and higher
        #     args.construct,
        #     genome,
        #     "-o", intra_sites  # available in minimap2 2.16 and higher
        #     ],
        #     check=True, capture_output=True
        # )
        minimap2_sp(args.construct, genome, intra_sites, args, secondary="no")
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "minimap2-search for identification of intra-sites failed"
            )
    out_size = os.path.getsize(intra_sites)
    if out_size == 0 and not checking:
        ilog.vlprint(
            "Search for intra-sites in genome has not yielded any results. If"
            " you are certain the reference contains preexisting integrations,"
            " Make sure to have the right integration sequence and to not use a"
            " genome already masked for your target. Otherwise this means the"
            " provided integration is just not present in the reference. In"
            " this case simply rerun intloc without the --intra flag.", 1
            )
        sys.exit()
    return intra_sites


def mask_genome(genome, intras_df, postfix=""):
    """Generate a genome with all integration sequences masked"""
    ilog.vlprint(
        f"Generating genome with masked intra-sites from {genome}",1
        )
    intras_dct = defaultdict(list)
    for index, intra in intras_df.iterrows():
        intras_dct[intra["qname"]].append((intra["qstart"], intra["qend"]))
    for chrom in intras_dct:
        intras_dct[chrom] = sorted(intras_dct[chrom], key=lambda x: x[0])

    genome_path, genome_name = os.path.split(genome)
    masked_name = ".".join(genome_name.split(".")[:-1]) + f"_{postfix}"
    masked_path = os.path.join(genome_path, masked_name)
    with open(genome, "r") as gen_in, open(masked_path, "w") as out:
        lc = 0
        curr_lines = list()
        contigs = dict()
        masked = 0
        for line in gen_in:
            lc += 1
            if line.startswith(">"):
                if lc > 1:
                    contigs[header].seq = curr_lines
                    curr_lines = list() 
                    masked += contigs[header].mask_seq(logger=ilog)
                    contigs[header].write_out(out)
                # minimap will only use the string before first " ", thus
                # the header must be reduced to this piece excluding ">"
                header = line.strip()
                mm_header = header.split(" ")[0][1:]
                seqs_to_mask = intras_dct[mm_header]
                tig = Contig(header, seqs_to_mask)
                contigs[header] = tig
            else:
                curr_lines.extend(list(line.strip()))
        contigs[header].seq = curr_lines
        masked += contigs[header].mask_seq(logger=ilog)
        contigs[header].write_out(out)
    ilog.vlprint(f"Masked {masked} sequences", 4)
    return masked_path


def mask_intra_sites(genome, args, chk=False, postfix="masked.tmp"):
    intra_sites = find_ref_intra_sites(genome, args, checking=chk)
    intras_df, cols = make_mm_paf_df(intra_sites)
    masked_genome = mask_genome(genome, intras_df, postfix=postfix)
    return masked_genome, intras_df, cols


def find_and_mask_all_intra_sites(args):
    """Find and mask preexisting intra-sites in genomic reference"""
    ilog.vlprint("Find and mask preexisting intra-sites in genomic reference",2)
    intra_name = os.path.split(args.construct)[1]
    mod_constr = "_".join(("_".join(intra_name.split(".")[:-1]).split(" ")))
    pfx = f'{mod_constr}_{args.ric}.fasta'
    genome_path, genome_name = os.path.split(args.genome)
    masked_name = ".".join(genome_name.split(".")[:-1]) + f"_masked_{pfx}"
    final_masked = os.path.join(genome_path, masked_name)
    intra_map_db = ".".join(final_masked.split(".")[:-1]) + ".idb"
    if os.path.exists(final_masked):
        ilog.vlprint(
            "Genome masked for intra-sites already exists. Skip masking step", 1
            )
        masked_gen2 = final_masked
        try:
            with open(intra_map_db, "r") as idb:
                data = idb.read()
                # len_filt_intras = json.loads(data)["len_filt_intras"]
                filt_dct = json.loads(data)["all_masked"]
                len_filt_intras = int_len_filter(filt_dct, args)
                ilog.vlprint(
                    f"Masked-intra-site-db (.idb) {masked_name} contains a tota"
                    f"l of {len(filt_dct)} records of masked sequences, with "
                    f"{len(len_filt_intras)} sequences of at least {args.ric}% "
                    f"of the length of the complete integration sequence.", 5
                    )
        except:
            ilog.vlprint(
                f"Masking step could not be skipped because {intra_map_db}"
                f" does not exist. Masking will be redone", 1
                )
            os.remove(final_masked)
            len_filt_intras, masked_gen2 = find_and_mask_all_intra_sites(args)
    else:
        # intra_sites = find_ref_intra_sites(args.genome, args)
        # intras_df, cols1 = make_mm_paf_df(intra_sites)
        # masked_genome1 = mask_genome(args.genome, intras_df, postfix="masked.tmp")
        masked_gen1, intras_df, cols1 = mask_intra_sites(args.genome, args)

        # postmask_intras = find_ref_intra_sites(masked_genome1, args, checking=True)
        # pm_intras_df, cols2 = make_mm_paf_df(postmask_intras)
        # masked_genome2 = mask_genome(masked_genome1, pm_intras_df, postfix=f"{pfx}")
        # os.remove(masked_genome1)
        masked_gen2, pm_intras_df, cols2 = mask_intra_sites(masked_gen1, args, chk=True, postfix=pfx)
        os.remove(masked_gen1)
        # masked_gen3, pm_intras_df3, cols3 = mask_intra_sites(masked_gen2, args, chk=True)
        
        # for gen in [masked_gen1]:
        #     os.remove(gen)

        filt_dct = filt_sort_paf_df(intras_df, cols1, no_culling=True, mm2_mapq=args.mm2_mapq)
        pm_filt = filt_sort_paf_df(pm_intras_df, cols2, no_culling=True, mm2_mapq=args.mm2_mapq)
        filt_dct.extend(pm_filt)
        # print("len(filt_dct)", len(filt_dct))
        len_filt_intras = int_len_filter(filt_dct, args)
        # print("len(len_filt_intras)", len(len_filt_intras))
        with open(intra_map_db,"w") as idb:
            # lfi_dct = {"len_filt_intras": len_filt_intras}
            # idb.write(json.dumps(lfi_dct, indent=4))
            all_masked = {"all_masked": filt_dct}
            idb.write(json.dumps(all_masked, indent=4))
        ilog.vlprint(
                    f"Masked-intra-site-db (.idb) {masked_name} contains a tota"
                    f"l of {len(filt_dct)} records of masked sequences, with "
                    f"{len(len_filt_intras)} sequences of at least {args.ric}% "
                    f"of the length of the complete integration sequence.", 5
                    )
        # ilog.vlprint(f"intras_df {len(intras_df)}", 2)
        ilog.vlprint(f"Masked {len(filt_dct)} sequences in reference genome", 1)
        ilog.vlprint(
            f"{len(len_filt_intras)} of the masked sequences have a length of"
            f" at least {args.ric}% of the integration sequence", 2
            )
        # ilog.vlprint(f"preex_intra_objs {len(len_filt_intras)}", 2)
    return len_filt_intras, masked_gen2, len(filt_dct)


def int_len_filter(filt_dcts, args):
    constr_len = int(filt_dcts[0]["tlen"]) if len(filt_dcts)>0 else 0
    args.constr_len = constr_len
    req_maplen = constr_len * (args.ric / 100)
    len_filt_intras = [i for i in filt_dcts if i["map_len"] >= req_maplen]
    return len_filt_intras


def gen_preex_intra_sites(site_dcts_lst, args):
    """Generate integration site objects for pre-existing intra-sites"""
    ilog.vlprint(
        f"Generating integration site objects for {len(site_dcts_lst)}"
        f" pre-existing intra-sites", 1
        )
    preex_intras = defaultdict(dict)
    for site in site_dcts_lst:
        new_intra = PotIntSite(
                                site["qname"],
                                (site["qstart"], site["qend"]),
                                [],
                                "n.a.",
                                args
                                )
        new_intra.preexisting_intra = True
        new_intra.preex_intra_alnm = site
        new_intra.approx_loc = new_intra.chrom_coord[0]
        new_intra.name = new_intra.chrom + "_" + str(new_intra.approx_loc) 
        preex_intras[new_intra.chrom][new_intra.approx_loc] = new_intra
    return preex_intras


def sort_cros_to_preex_intras(preex_intras, filt_rtr_als, cand_reads):
    """Sort candidate-read-objects first to read-to-ref alignments, then"""
    """ those in turn to preexisting-intra-site objects"""
    ilog.vlprint(
        "Sorting candidate read alignments to preexisting intra-sites", 1
        )

    amended_als = amend_mm_read2ref_als(filt_rtr_als, cand_reads)
    filt_rtr_als = {d["qname"]:d for d in amended_als}

    attach_stats = defaultdict(int)
    not_preexisting = []
    for al in filt_rtr_als:
        al = filt_rtr_als[al]
        tname, qname = al["tname"], al["qname"]
        attached_to = 0
        for prexi_pos in preex_intras[tname]:
            # if abs(prexi_pos - al["tstart"]) < al["qlen"]:
            if abs(prexi_pos - al["int_pos_chrom"]) < 10000:
                preex_intras[tname][prexi_pos].minimap_supp[qname] = al
                attached_to += 1
        if attached_to == 0:
            not_preexisting.append(al)
        attach_stats[attached_to] += 1
    ilog.vlprint(f"rtr-aln to preex-intra  attachment stats: {attach_stats}", 1)
    return not_preexisting


def run_il_intra(args):
    args.mm2_mapq = -1
    start = report_time("il_intra", end=False)
    len_filt_intras, masked_genome, msk_ct = find_and_mask_all_intra_sites(args)
    preex_intra_objs = gen_preex_intra_sites(len_filt_intras, args)
    cr_mapped = mm2_cand_reads(args.construct, args.reads_file, args.outdir, args)
    if cr_mapped.split(".")[-1]=="sam":
        cr_found_ct = extract_candidate_reads(cr_mapped)
        cr_paf = convert_sam2paf(cr_mapped)
    else:
        cr_paf = cr_mapped
        search_list = get_qnames_from_sam_or_paf(cr_mapped)
        crf_counter = get_reads_by_header(
            search_list, args.reads_file, "Integration_candidate_reads.fasta"
            )
        cr_found_ct = [k for k in crf_counter if crf_counter[k]>0]
    ilog.vlprint(f"Extracted {len(cr_found_ct)} integration candidate reads", 0)

    cand_df, cols = make_mm_paf_df(cr_paf, quick_mm2=True)
    cand_dcts = filt_sort_paf_df(cand_df, cols, mm2_mapq=args.mm2_mapq, no_culling=False)
    ilog.vlprint(f"Retained {len(cand_dcts)} reads after mapq filtering", 0)
    len_filt_crds = int_len_filter(cand_dcts, args) if len(cand_dcts)>0 else []
    ilog.vlprint(f"Retained {len(len_filt_crds)} reads after len filtering", 0)

    with open("candidate_reads_filt.tmp", "w") as tmp:
        for dct in len_filt_crds:
            tmp.write(f"{dct['qname']}\t{dct['qlen']}\n")
    cand_read_objs = gen_mmap_candread_inst(len_filt_crds, args)
    ilog.vlprint(f"Generated {len(cand_read_objs)} can_read objects", 0)
    write_reads_wo_integration(cand_read_objs)

    full_int_len, trunc_constr_len, prov_constr_len = write_intread_summary_txt(cand_read_objs, args)

    int_len = full_int_len if len(full_int_len)>0 else trunc_constr_len
    int_len = sum(int_len)/len(int_len) if len(int_len)>0 else 0
    if len(full_int_len) > 0:
        ilog.vlprint(
            f"Avg len of full integrations: "
            f"{sum(full_int_len)/len(full_int_len)}", 0
            )
    if len(trunc_constr_len) > 0:
        ilog.vlprint(
            f"Avg len of truncated integrations: "
            f"{sum(trunc_constr_len)/len(trunc_constr_len)}", 0
            )
    if all([len(full_int_len)==0, len(trunc_constr_len)==0]):
        ilog.vlprint(f"No integration candidate reads were found.", 0)

    read2ref_paf = mmquick_int_genome(masked_genome, args.outdir, args)
    filt_rtr_als = filt_mm_rtr_als(read2ref_paf, args)
    not_preex = sort_cros_to_preex_intras(preex_intra_objs, filt_rtr_als, cand_read_objs)
    pre_ios_lst = [pio for pos in preex_intra_objs.values() for pio in pos.values()]
    pre_ios_confirmed = len([pio for pio in pre_ios_lst if len(pio.minimap_supp) > 0])
    pre_ios_unsupported = len([pio for pio in pre_ios_lst if len(pio.minimap_supp) < 1])
    opt_info = {
        "masked sequences similar to int.-seq.": msk_ct,
        f"required int. cont. (ric) [percent of full length]": args.ric,
        f"Integrations (>=ric) in genomic reference": len(pre_ios_lst),
        "confirmed intragenic sites in reference": pre_ios_confirmed,
        "unsupported intragenic sites in reference": pre_ios_unsupported
    }
    preex_int_dict = write_int_csv(pre_ios_lst, prefix="preex_")
    draw_generic_ideogram(preex_int_dict, pre_ios_lst, args, prefix="preex_")
    if args.species:
        args.species = get_full_species_name(args)
        line_spacing = {
            "Homo sapiens": 5,
            "Saccharomyces cerevisiae": 1,
            "Mus musculus": 5,
            "Arabidopsis thaliana": 0,
            "Rattus norvegicus": 0,
            "Drosophila melanogaster": 0,
            "Escherichia coli": 0,
            "Danio rerio": 0,
            "Caenorhabditis elegans": 0
            }
        generate_spec_ideogram(line_spacing, args, prefix="preex_")
    pot_int_objs = mmquick_gen_int_sites(not_preex, args)
    ilog.vlprint(
        f"Generated {len(pot_int_objs)} potential novel integration sites", 1
        )
    reads_info = {}
    pot_int_objs = merge_pot_int_sites(pot_int_objs, n=1000, iteration=1)[0]
    ilog.vlprint(f"{len(pot_int_objs)} pot.-int.-sites remain after merging",1)
    pot_int_objs = filter_potintsites(pot_int_objs, args.min_supp, reads_info, args)
    opt_info["int sites detected by sequencing"] = len(pot_int_objs) + pre_ios_confirmed
    ilog.vlprint(f"{len(pot_int_objs)} PotIntSites remain after filtering",1)
    int_dict = write_int_csv(pot_int_objs)
    draw_generic_ideogram(int_dict, pot_int_objs, args)
    end = report_time("il_intra", end=True)
    intcnt = len(pot_int_objs)
    mmer_dct = contains_multimers(cand_read_objs)[1]
    av_crl = report_reads_info(reads_info, cr_found_ct, cand_read_objs, mmer_dct, full_int_len, trunc_constr_len, prov_constr_len, args)[1]
    multmers = None
    # Detection of multimers not attempted: the cand-read-method
    # check_for_concatemers() has not been tested with MmapCandidateRead yet,
    # although it has been inherited from BlastCandidateRead and theoretically 
    # could/should work; currently Evidence for multimers is not reported
    # for il_quick_mm2 and il_intra
    # multmers = "n.a."

    return pot_int_objs, av_crl, start, end, intcnt, multmers, reads_info, opt_info