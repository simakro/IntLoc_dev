import os
import pandas as pd
import subprocess
from collections import defaultdict

from il_subprocesses import caplog_sp_error, minimap2_sp
from il_classes import PotIntSite
from il_evaluate import make_mm_paf_df, filt_sort_paf_df
from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

def mmap_int_genome(genome, outdir, pot_int_sites, args):
    """Polyclonal mode: Identify potential integration locations by mapping  
    candidate reads against the reference using minimap2."""

    ilog.vlprint(
        'Polyclonal mode: Starting minimap2-search to identify location of '
        'integrations in reference genome', 1
        )
    os.chdir(outdir)
    if args.no_culling:
        int_cand_reads = "Integration_candidate_reads.fasta"
    else:
        int_cand_reads = "cand_reads_without_integration.fasta"
    
    outfile = os.path.join(outdir, "Integration_sites_minimap2")
    try:
        # subprocess.run(
        #         [
        #         args.minimap2,
        #         "-N", "0",
        #         "-x", f"map-{args.seq_tec}",
        #         genome,
        #         int_cand_reads,
        #         "-o", outfile
        #         ],
        #     check=True, capture_output=True
        #     )
        minimap2_sp(genome, int_cand_reads, outfile, args, secondary="no")
        with open(outfile, "r") as out:
            data = out.read().split("\n")
            refhit_mappings = [l.split("\t")[0] for l in data]
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
    
    try:
        mm_df, cols = make_mm_paf_df(outfile)
        mm_alnm = filt_sort_paf_df(mm_df, cols, mm2_mapq=0, no_culling=args.no_culling)
        rhr_after_filt = set([aln["qname"] for aln in mm_alnm])
        ilog.vlprint(
            f"{len(rhr_after_filt)} mm2-aligned reads remain after filtering",4
            )

        sort_count :dict= defaultdict(list)
        aln_ct = len(mm_alnm)
        absorbed = []
        
        for site in pot_int_sites:
            for aln in mm_alnm:
                if aln["qname"] in site.read_names:
                    site.minimap_supp[aln["qname"]] = aln
                    sort_count[aln["qname"]].append(site.name)
                    absorbed.append(aln["qname"])
            mm_alnm = [al for al in mm_alnm if al not in absorbed]
            
        overhang = {
            aln["qname"]:aln for aln in mm_alnm if aln["qname"] not in absorbed
            }
        ilog.vlprint(
            f"{len(set(overhang))} mm2-mapped reads have not been sorted to sites", 2
            )
        for aln in overhang:
            ilog.vlprint(aln, 6)
        sort_al = len(sort_count)
        ilog.vlprint(
            f"{sort_al}/{aln_ct} mm2-alignments were sorted to PotIntSite-objs",
            4
            )
        for al in sort_count:
            if len(sort_count[al]) > 1:
                ilog.vlprint(
                    f"{al} has been sorted to {len(sort_count[al])} different"
                    f" potential Integration sites. {sort_count[al]}", 2
                )
    except pd.errors.EmptyDataError as e:
        ilog.vlprint(
            "mmap_int_genome for identification of additional integration-sites"
            " in polyclonal mode failed", 0
            )

    return overhang.values()


def gen_mmonly_integration_sites(overhang, cand_read_objs, args):
    """Generate PotIntSite objects based on mm2 support only"""
    ilog.vlprint("Generate PotIntSite/s based on overhang mm2 alignment/s", 1)

    cros = {cro.RNAME:cro for cro in cand_read_objs}
    mm_only_pis = []
    for aln in overhang:
        piso = PotIntSite(
            aln["tname"], 
            (aln["tstart"], aln["tend"]),
            [],
            None,
            args
            )
        piso.minimap_supp[aln["qname"]] = aln
        piso.mm_only_cand_reads[aln["qname"]] = cros[aln["qname"]]
        piso.calc_mm2_loc()
        piso.re_calc_attributes()
        piso.rename()
        mm_only_pis.append(piso)
        ilog.vlprint(f"Generated mm-only PotIntSite {piso.name}", 3)
    return mm_only_pis


def run_il_polyclonal(genome, outdir, pot_int_sites, cand_read_objs, args):
    if not args.quick_mm2:
        overhang = mmap_int_genome(genome, outdir, pot_int_sites, args)
        mm_only_pis = gen_mmonly_integration_sites(overhang, cand_read_objs, args)
        pot_int_sites.extend(mm_only_pis)
    return pot_int_sites


        
    

