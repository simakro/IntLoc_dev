import os
import gzip
import json

import il_features as feat
import il_spec_ideo as ideo
from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

def write_reads_wo_integration(_cand_read_lst):
    """Excise integrations from cand-reads for mm2 evaluation (and ava-align)"""
    ilog.vlprint(
        "Excising integrations from cand-reads for evaluation, "
        "ava-align, quick_mm2, intra or polyclonal mode", 3
        )

    insuff_gen_seq = []
    with open("cand_reads_without_integration.fasta", "w") as f:
        for read in _cand_read_lst:
            if len(read.seq_wo_constr) > 18:
                f.write(">" + read.RNAME + "\n")
                f.write(read.seq_wo_constr.strip() + "\n")
            else:
                insuff_gen_seq.append(read)
    
    ilog.vlprint(
        f"{len(insuff_gen_seq)} reads not written to cand_reads_without_integra"
        f"tion.fa, for lack of genomic sequence after excission of int-seq.", 3
        )
    return insuff_gen_seq


def write_output_txt(pot_int_site_objs):
    """Generate int-report txt-file and export PotIntSite alignment-info as json"""

    with open('Integration_Report.txt', 'w') as rep, open('IntSitesAlignments.json', 'w') as alr:
        rep.write(
            f"A total number of {len(pot_int_site_objs)} potential integration "
            f"sites were identified based on mapping vs target genome.\n"
            )
        ilog.vlprint(
            f"A total number of {len(pot_int_site_objs)} potential integration sites "
            f"were identified based on mapping vs target genome.", 1
            )
        alignments_by_intsite = {}
        for int_site in pot_int_site_objs:
            rep.write(
                f"\tIntegration {int_site.name} is supported by {len(int_site.supp_als)} reads:\n"
                )
            ilog.vlprint(
                f"Integration {int_site.name} is supported by {len(int_site.supp_als)}"
                f" reads.",2
                )
            alignments_by_intsite.update(json.loads(int_site.to_json()))
            for al in int_site.supp_als:
                rep.write(f"\t\t {al.READ} {al.read_len}\n")
        alr.write(json.dumps(alignments_by_intsite, indent=2))
    return len(pot_int_site_objs)


def write_intread_summary_txt(gen_ca_re_in, args):
    int_constr_len = []  # lengths of full integrations found in reads
    trunc_constr_len = [] # lengths of truncated integrations found in reads
    mult_constr_len = [] # cumulated lengths of multiple int-seqs found in reads
    prov_constr_len = set()  # provided construct length
    with open("integration_read_summary.txt", "w") as sum_rep:
        for al in gen_ca_re_in:
            al.print_al_sum(sum_rep)
            int_compl = al.trunc_constr_read_pos()
            sum_rep.write(f"{al.RNAME}:{int_compl}\n")
            ilog.vlprint(f"{al.RNAME}:{int_compl}", 2, logging=False)
            prov_constr_len.add(al.constr_len)
            multi_read = any([al.multiple_sep_int, al.multimer_int])
            if not al.trunc_constr and not multi_read:
                int_constr_len.append(al.constr_len_in_read)
            elif multi_read:
                mult_constr_len.append(al.constr_len_in_read)
            else:
                trunc_constr_len.append(al.constr_len_in_read)
    if len(prov_constr_len)>0:
        prov_constr_len = sum(prov_constr_len)/len(prov_constr_len)
    else:
        with open(args.construct, "r") as cf:
            lines = [line.strip() for line in cf.read() if not line.startswith(">")]
            prov_constr_len = len("".join(lines))
    return int_constr_len, trunc_constr_len, prov_constr_len


def generate_spec_ideogram(line_spacing, args, prefix=""): # outdir, res_dir, species, 
    """Wrapper that integrates feature assignment and spec_ideogram creation"""
    ilog.vlprint("Initiating feature assignment and ideogram creation", 2)
    with open(os.path.join(args.res_dir, "genome_features.json"), "r") as fd:
            feat_dict = json.load(fd)
    features = gzip.open(os.path.join(args.res_dir, feat_dict[args.species]))
    # Integration_Report.csv is generated in il_report
    int_rep = f"{prefix}Integration_Report.csv"
    feat.il_features(
        features, os.path.join(args.outdir, int_rep), prefix=prefix
        )
    outfile = f"{prefix}ideogram.svg"
    try:
        ideo.run_spec_ideo(
            args.res_dir,
            args.outdir,
            outfile,
            line_spacing[args.species],
            args.species,
            args,
            prefix=prefix
            )
    except Exception as e:
        ilog.vlprint(
            f"WARNING: An unexpected error occured during construction"
            f" of species specific ideogram. {e}",
            0
        )

