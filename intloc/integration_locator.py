# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

from time import localtime, strftime
from collections import defaultdict, Counter
import subprocess
import argparse
import os
import platform
import sys
import gzip

sys.path.append(os.path.split(__file__)[0])
from version import __version__
from il_logging import Ilogger, initialize_loggers
from il_chk_conv import *
from il_report import *
from il_get import *
from il_asm_db_tools import *
from il_subprocesses import *
from il_pdf import write_pdf
from il_classes import *
from il_dumps import *
import il_mod_ideo
import il_html
import il_ava_clustering
import il_evaluate
import il_map_figs
import il_defaults


ilog = Ilogger()
ilog.module = __name__

DEFAULTS = il_defaults.DEFAULTS

def get_arguments():
    parser = argparse.ArgumentParser(
        description="Integration locator detects the number and location of "
        "integrations of transgenic constructs or other inserted DNA sequences "
        "within a genome using long-read sequencing data.",
        add_help=True
        )

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument(
        "-c", "--construct", required=True, 
        help="Path to fasta file containing the sequence of the integrating "
        "construct"
        )
    required_args.add_argument(
        "-o", "--outdir", required=True, help="Path for output directory"
        )
    
    read_input = required_args.add_mutually_exclusive_group(required=True)
    read_input.add_argument(
        "-f", "--reads_file",
        help="Path to a single fasta/q file containing all sequencing reads"
        )
    read_input.add_argument(
        "-d", "--reads_dir", help="Path to a directory containing multiple"
        " fasta/q files. Not recursive, subdirectories will be ignored."
        " (mutually exclusive with -f/--reads_file). Not compatible with"
        " --intra, --quick, --ava_clust and --polyclonal."
        )

    genome_choice = required_args.add_mutually_exclusive_group(required=True)
    genome_choice.add_argument(
        "-g", "--genome", help="Path to reference genome"
        )
    genome_choice.add_argument(
        "-w", "--dwnl_genome",
        help="Specify a species (full scientific name) for download of latest "
        "assembly (mutually exclusive with -g/--genome).",
        # a specific assembly can be specified using "Species name=asm_dir," 
        # w.o. "_genomic.fna.gz" postfix:
        # e.g. "Homo sapiens=GCF_009914755.1_T2T-CHM13v2.0"
        default=DEFAULTS["dwnl_genome"]
        )

    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument(
        "--spec", dest='species', default=DEFAULTS["species"], 
        help="Select species (e.g. \"human\") to enable output of species "
        "specific ideogram and annotations"
        )
    optional_args.add_argument(
        "--min_read_supp", dest="min_supp", default=DEFAULTS["min_supp"],
        help=f"Set minimum read support for integration validation"
        f" [{DEFAULTS['min_supp']}]", type=ranged_typechk(1,1000,int)
        )
    optional_args.add_argument(
        "--ava_clust", dest='cluster_reads', default=False, action='store_true', 
        help="Cluster reads based on all-vs-all alignment. Comparison of read-"
        "clusters to integration-sites identified by alignment of reads to the "
        "reference genome may assist identification of potential duplicates "
        "that could arise from integration into highly repetetive genomic loci."
        )
    optional_args.add_argument(
        "--req_int_cont", dest='ric', type=ranged_typechk(1,100,int),
        default=DEFAULTS["req_int_cont"],
        help="Only include sequencing reads containing at least x%% of the int"
        "egrating sequence as candidates for the detection of integration sites"
        ". [default=20]"
        )
    # optional_args.add_argument(
    #     "--req_ref_cont", dest='rrc', default=20,type=ranged_typechk(1,100,int),
    #     help="Only include sequencing reads containing at least x%% of the int"
    #     "egrating sequence as candidates for the detection of integration sites"
    #     ". [default=20]"
    #     ) => it may still be possible to use this parameter using absolute 
    #   base-count of ref-portion instead of a proportion
    optional_args.add_argument(
        "--cov", dest='coverage', default=None,
        type=ranged_typechk(1,100000,float),
        help="Experimental. For whole genome sequencing data an approximate"
        " coverage value may be supplied as integer or floating point number"
        " through this argument. This will be used to calculate the minimum"
        " read support for integration validation."
        " Assumes relatively homogenous coverage (+/- 50%%) and diploidy."
        " Only recommended for cov > 15x. For data with uncertain or highly in-"
        " homogenous coverage, such as data generated using target-enrichment"
        " strategies, adaptive sampling or amplicon-sequencing, as well as for"
        " very low coverage WGS data, use the --min_read_supp arg to adjust"
        " the minimum required read support to an appropriate value manually."
        " [default=None]"
        )
    optional_args.add_argument(
        "--keep_nonint", dest='keep_nonint', default=False, action="store_true",
        help="Write non-int reads at int-sites (typically stemming from homolog"
        " without integration) to a fasta file for scrutiny [default=False]"
        )
    optional_args.add_argument(
        "--verbosity", dest='verbosity', default=1,
        type=ranged_typechk(-1,100,int), 
        help="Adjust level of output to the terminal"
        )
    optional_args.add_argument(
        "-v", "--version", action="version",
        version=f"integration_locator {__version__}",
        help="Print version and exit"
        )
    optional_args.add_argument(
        "--file_check", dest='check_reads', default="quick",
        help="Adjust thoroughness of file-check; quick/extensive [quick]"
        )
    optional_args.add_argument(
        "--ntigs", dest='ntigs', default=False, type=ranged_typechk(0,1000,int),
        help="Set number of contigs to be included in generic ideogram. If not"
        " set this number is handled by the Nx parameter. [Default=False]"
        )
    optional_args.add_argument(
        "--include_Nx", dest='Nx', default=98, type=ranged_typechk(1,100,float),
        help="Choose cutoff for contigs to be included in figures in gneric "
        "ideogram. Contigs comprising at least x percent of bases in RefGenome "
        "[98] will be included. This parameter can be helpful to avoid graphs "
        "becoming overcrowded when working with reference genomes that contain "
        "large amounts of alternative loci or a high number of miniscule "
        "unplaced contigs. Don't worry, this applies only to the figures 1-3 in"
        " the pdf-report. Thus, no sequences will be ignored during integration"
        " search, i.e. there is no risk of missing out on integration sites. "
        "Species specific ideogram is not affected."
        )
    optional_args.add_argument(
        "--bait", dest='bait', default=False,
        help="This argument can be supplied in addition to the construct/-c "
        "argument. It is required if the integration construct shares partial "
        "similarity with the reference sequence when intloc is not run in intra"
        "mode (--intra argument). Intloc will inform you if the use of a bait"
        " sequence is necessary. Provide the path to a fasta file "
        "containing only part of the integrating construct as bait for "
        "identification of integration evidence reads. The bait sequence must "
        "not share similarity with reference. [False]"
        )
    optional_args.add_argument(
        "--seq_tec", dest='seq_tec', default="ont",
        help="Specify (long-read) sequencing technology (ont/pb). "
        "[default='ont']"
        )
    optional_args.add_argument(
        "--min_dist_ints", dest='min_dist_ints', default=25,
        type=ranged_typechk(10,50000,int),
        help="Sets the minimimun distance in bp that is required in order to "
        " count integrations on the same DNA molecule as separate entities or"
        " as part of a single (potentially complex) integration event."
        " [defaults by mode: standard=25, quick_mm2=1000]"
        )
    optional_args.add_argument(
        "--untrimmed", dest='untrimmed', default=0,
        type=ranged_typechk(0,250,int),
        help="If you are using untrimmed raw-reads as input, you can provide "
        " the expected lengths of sequencing adapters + barcodes (if used) at"
        " the ends of each read via this parameter. [default=0]"
        )
    optional_args.add_argument(
        "--polyclonal", dest='polyclonal', default=False, action='store_true',
        help="Use this flag to set min_read_supp=1, prevent merging of int-sites"
        " that are spaced less than min_dist_ints from each other and prevent"
        " exclusion of sites that share extensive similarity and overlap."
        " [default=False]"
        )
    optional_args.add_argument(
        "--intra", dest='intra', default=False, action='store_true',
        help="Intra-mode allows to find pre-existing as well as newly occuring"
        " integrations of a provided query sequence in the genome. For example"
        " this mode can be used to identify the locations of LINEs in a genome."
        " [default=False]"
        )
    optional_args.add_argument(
        "--map_qual", dest='map_qual', default=90,
        help="Specify cutoff for required completeness (not identity) of read"
        " mappings to the reference. If the genome of your sample does not"
        " deviate significantly from the reference, this value should be"
        " 3-10 %% lower than the accuracy of the sequencing reads."
        " Increasing the tolerance by lowering this value may help to reveal" 
        " integrations in genomic locations of your sample that significantly"
        " deviate from the reference. However, it may increase the risk"
        " for false positives as well. [default=90]"
        )
    optional_args.add_argument(
        "--no_multi_ints_info", dest='no_multi_ints', default=False,
        action='store_true',
        help="Toggle precise analysis of reads with multiple insertions of the "
        "integrating sequence. When deactivated those reads will typically just"
        " support one integration site and linkage analysis does not take place"
        ". [default=False]"
        )
    optional_args.add_argument(
        "--cleanup", dest='cleanup_aux', default=True, action='store_false',
        help="Deactivate cleanup of temporary and auxiliary files."
        " Relevant for debugging only."
        )
    optional_args.add_argument(
        "--skip_eval", dest='skip_eval', default=False, action='store_true',
        help="Skip evaluation module. Evaluation for the identification of false"
        " positives can significantly increase the runtime (in dependence of"
        " genome size and number of integrations). You can skip it when in a"
        " hurry and do not fear no FPs."
        )
    optional_args.add_argument(
        "--quick", dest='quick_mm2', default=False, action='store_true',
        help="Quick and dirty integration detection (using only minimap2) witho"
        "ut thorough evaluation of sites. Often sufficient for easy cases (smal"
        "ler non-repetetive genomes, low number of integrations), but typically"
        " produces poorer results in more demanding cases (large complex genome"
        "s, high number of integrations) than standard mode."
        "Not compatible with following arguments: --reads_dir, --ava_clust"
        )
    optional_args.add_argument(
        "--int_len", dest='int_len', default=False,
        type=ranged_typechk(25,10000000,int),
        help="Manually set the expected length of integrations. For WGS data"
        " this is not required nor recommended, since IntLoc automatically"
        " determines the lengths of integrations in the genome. However, this"
        " parameter can be helpful when only using part of the integration as"
        " bait, or when using enrichment methods (e.g. Cas9) that always result"
        " in truncated integration sequences in the reads. [default=False]"
        )
    optional_args.add_argument(
        "--cores", dest='cores', default=False,
        type=ranged_typechk(1, 1024, int, out_type=str),
        help="Specify number of cores to use. Intloc by default attempts to"
        " run parallelized tasks on all but two of the cores available on the"   
        " current machine (e.g. using 6 of 8 available cpu cores)."
        )
    optional_args.add_argument(
        "--gene_dens", dest='gene_dens', default=False, action='store_true',
        help="Derive ideograms with gene-density information."
        )
    optional_args.add_argument(
        "--feat_dens", dest='feat_dens', default=False,
        help="Derive ideograms with feature-density information. A feature clas"
        "s and subclass have to be provided separated by a comma. (--feat_dens "
        "gene,protein_coding reproduces the function of the --gene_dens flag)"
        )
    optional_args.add_argument(
        "--color_ints", dest='color_ints', default="lime",
        help="Specify color of int labels in figures. Possible values: gre" 
        "en, blue, cyan, magenta, yellow, lime, olive, orange, pink, purple,"
        " teal, violet [default='lime']"
        )
    optional_args.add_argument(
        "--override", dest='override', default=False, action='store_true',
        help=argparse.SUPPRESS
        )
    optional_args.add_argument(
        "--no_culling", dest='no_culling', default=False, action='store_true',
        help=argparse.SUPPRESS
        )
    optional_args.add_argument(
        "--gui", dest='gui', default=False, action='store_true', 
        help=argparse.SUPPRESS
        )
    optional_args.add_argument(
        "--testing", dest='testing', default=False, 
        # help="Provide a csv containing information about known integration"
        # " sites for comparison to sites identified by IntLoc. csv-format:" 
        # " ID,chromosome,bp-location. No header.", 
        help=argparse.SUPPRESS
        )
    # optional_args.add_argument(
    #     "--test_batch", dest='test_batch', default=False, 
    #     # help="Run tests for a number of commands and collect results." 
    #     help=argparse.SUPPRESS
    #     )
    optional_args.add_argument(
        "--cat_log", dest='cat_log', default=False, action="store_true",
        # help="Print intlog.log to stdout"
        help=argparse.SUPPRESS
        )
    optional_args.add_argument(
        "--dep_ver", dest='dep_ver', default={},
        # help="dict for storage of infos about dependency versions"
        help=argparse.SUPPRESS
        )
    # optional_args.add_argument(
    #     "--env", dest='env', default=None,
    #     # help="dict to be passed to subprocesses env kwarg in intloc_self_sp"
    #     help=argparse.SUPPRESS
    #     )

    args = parser.parse_args()

    return args   


def gen_output_dir(outdir, reads_dir=False, override=False, gui=False):
    if os.path.exists(outdir) and not reads_dir:
        ilog.vlprint('Output directory already exists.', 10)
        if gui:
            answer = "N"
        else:
            answer = input(
                'Output directory already exists. Do you want to proceed? (Y/N)'
                ).upper()
        if answer == 'N':
            ilog.vlprint(
                'Declined writing into existing output directory. Aborting', 10
                )
            sys.exit()
        elif answer == 'Y':
            ilog.vlprint('Accepted writing into existing output directory.', 10)
            pass
        else:
            ilog.vlprint('Invalid entry. Aborting program', 0)
            sys.exit()
    else:
        try:
            if override:
                os.makedirs(outdir, mode=777, exist_ok=True)
                os.chmod(outdir, 0o777)
            else:
                os.makedirs(outdir, mode=777, exist_ok=False)
                os.chmod(outdir, 0o777)
        except OSError:
            ilog.vlprint(
                "Output directory can't be created. Either path does not exist,"
                " or output directory already exists. Please choose a different"
                " path/name for output directory or rename/delete existing"
                " directory.", 0
                )
            sys.exit()
        ilog.vlprint('Created output directory.', 10)


def check_constr_for_refseqs(genome, construct, args, bait=False):
    """Check if the integrating construct contains sequences with significant
        similarity to sequences in the reference genome."""

    def search_non_ref_seq(const_chk_file):
        with open(const_chk_file, "r") as in_file:
            ref_ranges_in_constr = []
            for line in in_file:
                if line.startswith("RNAME"):
                    pass
                else:
                    ls = line.split(",")
                    ref_range = (int(ls[5]), int(ls[6]))
                    rr_sort = tuple(sorted(ref_range))
                    ref_ranges_in_constr.append(rr_sort)
            ref_ranges_sorted = sorted(ref_ranges_in_constr, key=lambda x: x[0])
            if len(ref_ranges_sorted)==1:
                aln = ref_ranges_sorted[0]
                if (aln[1]-aln[0]+1) == int(ls[1]):
                    ilog.vlprint(
                        "The construct sequence appears to entirely consist of"
                        " reference sequence. Try using --intra.", 1
                        )
                    sys.exit()
            else:
                gaps = []
                for i in range(len(ref_ranges_sorted) - 1):
                    dist = abs(ref_ranges_sorted[i][1] - (ref_ranges_sorted[i+1][0] - 1))
                    gaps.append(tuple([i, dist]))
                max_gap = max(gaps, key=lambda x: x[1])
                sugg_range_start = ref_ranges_sorted[max_gap[0]][1] + 1
                sugg_range_stop = ref_ranges_sorted[max_gap[0] + 1][0] - 1
        return sugg_range_start, sugg_range_stop


    if bait:
        construct = bait
        ilog.vlprint(
            'Checking if the bait contains sequences with significant similarity' 
            '  to sequences in the reference genome.', 1
            )
    else:
        ilog.vlprint(
            'Checking if the integrating construct contains sequences with ' 
            'significant similarity to sequences in the reference genome.', 1
            )

    os.chdir(args.outdir)
    genome = check_for_refdb(genome, args)
    constr_chk = os.path.join(args.outdir, "constr_chk_refseq.out")

    tmp = "constr_chk.tmp"
    outfile = open(constr_chk, 'a')
    outfile.write(
        'RNAME, qlen, CHR, slen, length, qstart, qend, sstart, send, pident,'
        ' mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus\n'
        )
    try:
        subprocess.run(
            [
            args.blastn,
            "-db", genome,
            "-query", construct,
            "-out", tmp,
            "-outfmt", "10 qaccver qlen saccver slen length qstart qend sstart"
            " send pident mismatch gapopen evalue bitscore qcovhsp qcovus",
            "-num_threads", args.cores,
            "-culling_limit", "1",
            "-evalue", "0.0001",
            "-perc_identity", "95"
            ],
            check=True, capture_output=True
            )
    except subprocess.CalledProcessError as e:
        caplog_sp_error(
            e, "BLAST search for reference sequences in integration construct"
            " failed"
            )

    with open(tmp, 'r') as temp:
        tmp_out = temp.read()
        outfile.write(tmp_out)
    outfile.close()
    os.remove(tmp)
    with open(constr_chk, 'r') as results:
        matches = []
        for line in results:
            if line.startswith("RNAME"):
                pass
            else:
                ls = line.split(",")
                matches.append(int(ls[4]))
        if len(matches) > 0:
            if sum(matches) > 18:
                sugg_range = search_non_ref_seq(constr_chk)
                ilog.vlprint(
                    f'The search sequence (integration construct/bait) contains'
                    f' sequences with significant similarity to native sequences' 
                    f' in the reference genome. This can negatively affect the '
                    f'accurate detection, quantification and placement of '
                    f'integrations. If you want to detect integrations (old and' 
                    f' new) of sequences known to be already present in the'
                    f' reference genome, please use the --intra argument.'
                    f' Otherwise, you can opt to use parts of the integration-'
                    f'construct, which do not occur in the reference, as bait. '
                    f'This will help to avoid false positive read hits and '
                    f'enable accurate integration identification.\n'
                    f'#####-IntLoc suggests to use bases {sugg_range[0]} to '
                    f'{sugg_range[1]} of the integration construct as bait-'
                    f'sequence. Place this sequence in a separate fasta file and'
                    f' add the --bait flag with the path to this file to your ' 
                    f'command.-#####', 1
                    )
                sys.exit()
        else:
            ilog.vlprint(
                'The search sequence (integration construct/bait) does not share'
                ' significant similarity to native sequences in the reference '
                'genome.', 2
                )

    return args.outdir, constr_chk, genome


def blast_cand_reads(construct, reads, args, reads_dir=False, bait=False):
    """Identify candidate reads containing integration of interest in sequencing
     data"""

    readpath = os.path.split(reads)
    workpath = readpath[0]
    os.chmod(workpath, 0o777)

    db_ext = ['.ndb', '.nal', '.nsq', '.00.nsq']
    if any([os.path.exists(reads + ext) for ext in db_ext]):
        ilog.vlprint('Read-db already exists. Skipping generation of read-db',1)
    else:
        ilog.vlprint('Generating BLAST-DB from reads file', 0)
        try:
            subprocess.run(
                [
                args.makeblastdb,
                "-in", reads,
                "-parse_seqids",
                "-dbtype", "nucl"
                ],
                check=True, capture_output=True
                )
        except subprocess.CalledProcessError as e:
            caplog_sp_error(
                e, "makeblastdb for identification of candidate reads failed"
                )
    os.chdir(args.outdir)

    if bait:
        cand_reads = os.path.join(args.outdir, "bait_candidate_reads")
    else:
        cand_reads = os.path.join(args.outdir, "candidate_reads")

    if os.path.exists(cand_reads) and not reads_dir:
        ilog.vlprint(
            'candidate_reads file exists. '
            'Skipping search for reads containing construct', 1
            )
    else:    
        tmp = "out.tmp"

        outfile = open(cand_reads, 'a')
        outfile.write(
            f'RNAME, qlen, CHR, slen, length, qstart, qend, sstart, send,  pident,'
            f' mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus, sstrand\n'
            )
        ilog.vlprint(
            'Starting BLAST-search to identify reads containing '
            'integration-construct sequences', 1
            )
        try:
            subprocess.run([
                args.blastn,
                "-db", reads,
                "-query", construct,
                "-out", tmp,
                "-outfmt", 
                "10 qaccver qlen saccver slen length qstart qend sstart send "
                "pident mismatch gapopen evalue bitscore qcovhsp qcovus sstrand",
                "-num_threads", args.cores,
                "-evalue", "0.000000000000000001",
                "-perc_identity", "80",
                "-dust", "no",
                "-max_target_seqs", "10000000"],
                check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            caplog_sp_error(
                e, "BLAST-search for identification of candidate reads failed"
                )

        with open(tmp, 'r') as temp:
            tmp_out = temp.read()
            outfile.write(tmp_out)
        outfile.close()

    return args.outdir, cand_reads


def blast_output_filter(input_blast_report, reads_dir=False, req_int_cont=20):
    """Filter candidate reads for results that contain at least {req_int_cont}%  
    of the integration of interest"""
    ilog.vlprint("Filtering BLAST-results", 4)

    ricf = 100/req_int_cont
    if os.path.exists(input_blast_report + '.fbr') and not reads_dir:
        filtered_blast_report = input_blast_report + '.fbr'
    else:
        filtered_blast_report = input_blast_report + '.fbr'
        fbr_out = open(filtered_blast_report, "w")

        with open(input_blast_report) as btr:
            for line in btr:
                if line.startswith('RNAME'):
                    fbr_out.write(line)
                else:
                    l = line.split(',')
                    if int(l[4]) > int(l[1])/ricf:
                        fbr_out.write(','.join(l))
                    else:
                        continue
        
        fbr_out.close()
    return filtered_blast_report


def retrieve_cand_reads(blast_rep_all, blast_rep_filt, reads, prefix=""):
    """Retrieve candidate reads containing integration from sequencing data"""
    ilog.vlprint('Retrieve candidate reads', 1)

    def retrieve_reads(blast_report, prefix=""):
        search_list = set()
        with open(blast_report, 'r') as fbr:
            for line in fbr:
                if line.startswith('RNAME'):
                    pass
                else:
                    l = line.split(',')
                    header = ">" + l[2]
                    search_list.add(header)
        
            cand_read_path = os.path.join(
                os.getcwd(), f"{prefix}_candidate_reads.fasta"
                )
        return search_list, cand_read_path

    sl_all, cr_path_all = retrieve_reads(blast_rep_all, prefix="Unfiltered")
    sl_filt, cr_path_filt = retrieve_reads(blast_rep_filt, prefix="Integration")
    
    all_count = get_reads_by_header(sl_all, reads, cr_path_all, specifier="all")
    filt_count = get_reads_by_header(sl_filt, reads, cr_path_filt, specifier="filtered")

    return cr_path_filt, filt_count, all_count


def gen_cand_read_instances(filtered_blast_report, args):
    """Generate candidate read objects"""
    ilog.vlprint("Generating candidate read objects", 2)

    cands = {}
    with open(filtered_blast_report) as fbr:
        for line in fbr:
            if line.startswith('RNAME'):
                pass
            else:
                l = line.strip().split(',')
                l = [
                    l[0],
                    int(l[1]),
                    l[2],
                    int(l[3]),
                    int(l[4]),
                    int(l[5]),
                    int(l[6]),
                    int(l[7]),
                    int(l[8]),
                    float(l[9]),
                    int(l[10]),
                    int(l[11]),
                    float(l[12]),
                    int(l[13]),
                    int(l[14]),
                    int(l[15]),
                    l[16]
                    ]
                ilog.vlprint(l, 3, logging=False)
                rname = l[2]
                if rname not in cands:
                    alignment = {
                                "rname": l[2],
                                "slen": l[3],
                                "length": l[4],
                                "crange": [l[5:7]],
                                "rrange": [l[7:9]],
                                "pident": l[9],
                                "mismatch": l[10],
                                "gapopen": l[11],
                                "evalue": l[12],
                                "bitscore": l[13],
                                "qcovhsp": l[14],
                                "qcovus": l[15],
                                "strand": {tuple(sorted(l[7:9])):l[16]},
                                "qlen": l[1]
                            }
                    cands[rname] = alignment
                else:
                    # if multiple alignments exist for one read add to or 
                    # average the respective values with the primary alignment
                    cands[rname]["length"] = cands[rname]["length"] + l[4]
                    cands[rname]["crange"].append(l[5:7])
                    cands[rname]["rrange"].append(l[7:9])
                    cands[rname]["pident"] = (cands[rname]["pident"] + l[9])/2
                    cands[rname]["mismatch"] = cands[rname]["mismatch"] + l[10]
                    cands[rname]["gapopen"] = cands[rname]["gapopen"] + l[11]
                    cands[rname]["evalue"] = (cands[rname]["evalue"] + l[12])/2
                    cands[rname]["bitscore"] = cands[rname]["bitscore"] + l[13]
                    cands[rname]["qcovhsp"] = cands[rname]["qcovhsp"] + l[14]
                    cands[rname]["qcovus"] = (cands[rname]["qcovus"] + l[15])/2
                    cands[rname]["strand"][tuple(sorted(l[7:9]))] = l[16]
    ilog.vlprint(f"\n{cands}\n", 4, logging=False)

    # First sort order of sublists in crange and rrange lists
    for k in cands:
        cands[k]["crange"] = sorted(cands[k]["crange"])
        cands[k]["rrange"] = sorted(cands[k]["rrange"])
    # Then sort values in those sublists
    for k in cands:
        n = 0
        for _ in cands[k]["crange"]:
            cands[k]["crange"][n] = sorted(cands[k]["crange"][n])
            n += 1
    for k in cands:
        n = 0
        for _ in cands[k]["rrange"]:
            cands[k]["rrange"][n] = sorted(cands[k]["rrange"][n])
            n += 1
    ilog.vlprint(f"{cands}\n", 3, logging=False)

    # Reduce alignment length by bases contained in overlaps or gaps between 
    # subalignments
    for k in cands:
        if len(cands[k]["crange"]) > 1:
            r = range(len(cands[k]["crange"])-1)
            for n in r:
                s = abs(cands[k]["crange"][n][1] - cands[k]["crange"][n+1][0])
                cands[k]["length"] = cands[k]["length"] - s
    ilog.vlprint(f"{cands}\n", 4, logging=False)        

    # For multiple subalignments fuse the cornering values of crange and 
    # rrange to show covered area
    for k in cands:
        # add orig_alignments to dct before fusing cornering vals
        cands[k]["orig_alnm"] = cands[k]["rrange"]
        if len(cands[k]["crange"]) > 1:
            r = len(cands[k]["crange"])-1
            cands[k]["crange"] = [cands[k]["crange"][0][0], cands[k]["crange"][r][1]]
            cands[k]["rrange"] = [cands[k]["rrange"][0][0], cands[k]["rrange"][r][1]]
        else:
            cands[k]["crange"] = cands[k]["crange"][0]
            cands[k]["rrange"] = cands[k]["rrange"][0]
    ilog.vlprint(f"{cands}\n", 4, logging=False)

    # Finally generate object instances of BlastCandidateRead-Class
    cand_read_lst = []
    ct = 0
    int_cand_reads = sl_fasta2dict("Integration_candidate_reads.fasta")
    for key in cands:
        ct += 1
        seq = int_cand_reads[f">{key}"]
        attrs = cands[key]
        inst_obj = BlastCandidateRead(
                                      ct, 
                                      attrs["rname"],
                                      attrs["slen"],
                                      attrs["length"],
                                      attrs["crange"],
                                      attrs["rrange"],
                                      attrs["pident"],
                                      attrs["mismatch"],
                                      attrs["gapopen"],
                                      attrs["evalue"],
                                      attrs["bitscore"],
                                      attrs["qcovhsp"],
                                      attrs["qcovus"],
                                      attrs["qlen"],
                                      seq,
                                      attrs["orig_alnm"],
                                      attrs["strand"],
                                      args
                                      )
        cand_read_lst.append(inst_obj)

    ilog.vlprint(f"Generated {len(cand_read_lst)} candidate read objects", 2)
    return cand_read_lst


def contains_multimers(_gen_ca_re_in):
    """Function that evaluates if multimers/concatmers of the inserting sequence
     have integrated into the genome, or if closely spaced integrations may have
     taken place."""
    ilog.vlprint("Checking for multimers and closely spaced integrations.", 2)

    multiple_ind = {}
    multimeric = {}
    simple = 0
    for cri in _gen_ca_re_in:
        if not any([cri.multiple_sep_int, cri.multimer_int]):
            simple += 1
        if cri.multiple_sep_int:
            multiple_ind[cri.RNAME] = cri
        elif cri.multimer_int:
            multimeric[cri.RNAME] = cri
        else:
            pass

    with \
        open("multimeric_int.fasta", "w") as mm_fa, \
        open("multiple_indep_ints.fasta", "w") as mi_fa, \
        open("Multimer_info.txt", "w") as mm_info:
        for read_name in multimeric:
            mm_fa.write(f">{read_name}\n")
            mm_fa.write(f"{multimeric[read_name].seq}\n")
            mm_info.write(
                f"Read {read_name} appears to contain multimeric insertions of "
                "the integrating DNA sequence\n"
                )
        for read_name in multiple_ind:
            mi_fa.write(f">{read_name}\n")
            mi_fa.write(f"{multiple_ind[read_name].seq}\n")
            mm_info.write(
                f"{read_name} appears to contain multiple individual insertions"
                f" of the integrating DNA sequence. Mutliple independent, but"
                f" very closely spaced insertions may be merged into a single"
                f" integration site, depending on the value of the min_dist_int"
                f" parameter.\n"
                )

    if any(multiple_ind) or any(multimeric):
        return True,  {"mult_ind": multiple_ind, "multimer": multimeric, "simple": simple}
    else:
        return False, {"mult_ind": {}, "multimer": {}, "simple": simple}


def blast_int_genome(genome, args):
    """Identify potential locations in the genome for the integrated sequences,
    by blasting candidate reads against the reference."""

    ilog.vlprint(
        'Starting BLAST-search to identify location of integrations in '
        'reference genome', 1
        )
    os.chdir(args.outdir)
    int_cand_reads = "Integration_candidate_reads.fasta"
    # Reads with excised integration are not intended to be used here, but in
    # il_ava_clustering. If used here int_site placement accuracy is reduced.
    # int_cand_reads = "cand_reads_without_integration.fasta"
    genome = check_for_refdb(genome, args)
    cand_ints = os.path.join(args.outdir, "Integration_sites")

    if os.path.exists(cand_ints):
        ilog.vlprint('File Integration_sites already exists.', 1)
        ilog.vlprint('Skipping search for candidate integrations in genome', 1)
    else:    
        tmp = "ints_out.tmp"
        outfile = open(cand_ints, 'a')
        outfile.write(
            'RNAME, qlen, CHR, slen, length, qstart, qend, sstart, send,'
            ' pident, mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus,'
            ' score, gaps, nident, positive, sstrand\n'
            )
        try:
            subprocess.run([
                args.blastn,
                "-db", genome,
                "-query", int_cand_reads,
                "-out", tmp,
                "-outfmt", "10 qaccver qlen saccver slen length "
                "qstart qend sstart send pident mismatch gapopen "
                "evalue bitscore qcovhsp qcovus score gaps nident "
                "positive sstrand", 
                "-num_threads", args.cores,
                "-culling_limit", "1",
                "-evalue", "0.0001",
                "-perc_identity", "80",
                "-gapopen", "0",
                "-gapextend", "0"],
                check=True, capture_output=True
                )
        except subprocess.CalledProcessError as e:
            caplog_sp_error(
                e, "BLAST search for identification of integration positions in"
                " reference failed"
                )

        with open(tmp, 'r') as temp:
            tmp_out = temp.read()
            outfile.write(tmp_out)
        outfile.close()

    return args.outdir, cand_ints, genome


def gen_read_to_ref_al(cand_ints, cand_read_lst, args):
    """LocalAlnm Version; Creation of read to reference alignment objects"""
    ilog.vlprint(f"Generating read-to-reference alignments", 3)
    al_dict = {}
    loc_al_ct = 0
    loc_al_reads = set()
    with open(cand_ints) as cin:
        for line in cin:
            if line.startswith('RNAME'):
                pass
            else:
                l = line.strip().split(',')
                ilog.vlprint(l, 3, logging=False)
                loc_al = LocalAlnm(*l)
                loc_al_ct += 1
                loc_al_reads.add(l[0])
                read = (loc_al.read, loc_al.read_len)
                chrom = (loc_al.chrom, loc_al.chrom_len)
                if read not in al_dict:
                    al_dict[read] = {chrom: [loc_al]}
                elif chrom in al_dict[read]:
                    al_dict[read][chrom].append(loc_al)
                else:
                    al_dict[read][chrom] = [loc_al]

    ilog.vlprint(f"al_dict: {al_dict}", 4, logging=False)
    ilog.vlprint(
        f"Found {loc_al_ct} local reference alignments"
        f" for {len(loc_al_reads)} reads", 8
        )

    # Generate object instances of ReadToRefAl-Class
    rtr_lst = []
    for read_key in al_dict:
        for read_obj in cand_read_lst:
            if read_obj.RNAME == read_key[0]:
                cand_read_obj = read_obj
                break
        try:
            inst_obj = ReadToRefAl(
                                    read_key[0],
                                    read_key[1],
                                    al_dict[read_key],
                                    cand_read_obj,
                                    args
                                    ) 
            rtr_lst.append(inst_obj)
            cand_read_obj.rtr = inst_obj
        except NameError as e:
            ilog.vlprint(str(e), 0)
            ilog.vlprint(
                f"In gen_read_to_ref_al cand_read_obj for {read_key[0]} "
                f"could not be found", 0
                )

    ilog.vlprint(f"rtr_lst: {rtr_lst}", 2, logging=False)
    ilog.vlprint(f"Generated {len(rtr_lst)} read-to-reference alignments", 3)
    return al_dict, rtr_lst


def nucleate_int_sites(rtr_lst, args):
    """Determine Integration sites from read to reference alignments"""
    ilog.vlprint(
        "Clustering read to reference alignments integration site nuclei ", 3
        )
    chr_dict = defaultdict(list)
    loc_dict = {}
    for al in rtr_lst:  
        if al.CHR in chr_dict:
            chr_dict[al.CHR].append(al)
        else:
            chr_dict[al.CHR] = [al]

    ilog.vlprint(f"chr_dict: {chr_dict}", 2, logging=False)
    for chrom in chr_dict:
        chr_dict[chrom] = sorted(chr_dict[chrom], key=lambda x: x.coord_evid)
        loc_dict[chrom] = {}
        for rtr_alnm in chr_dict[chrom]:
            inserted = False
            if len(loc_dict[chrom]) > 0:
                for nucl_site in loc_dict[chrom]:
                    inserted = loc_dict[chrom][nucl_site].compare_incoming_rtr(rtr_alnm)
                    if inserted:
                        break
            if not inserted:
                new_site = NucleationSite(chrom, rtr_alnm, args)
                loc_dict[chrom][new_site.id] = new_site
    
    nucl_ct = 0
    for chrom in loc_dict:
        for nucl_site in loc_dict[chrom]:
            nucl_ct += 1
            loc_dict[chrom][nucl_site].calc_metrics()
            loc_dict[chrom][nucl_site].report_nucl_site()

    ilog.vlprint(f'Generated {nucl_ct} nucleation sites.', 2)
    ilog.vlprint(f'loc_dict: {loc_dict}', 2, logging=False)
    return loc_dict


def gen_integration_sites(loc_dict, outdir, genome, args):
    """Determine Integration sites from read to reference alignments"""
    ilog.vlprint(
        "Generating Integration sites based on int-site nuclei", 3
    )
    # get blast parameters for candidate reads vs ref_genome
    try:
        chrom = list(loc_dict.keys())[0]
        sample_read = list(loc_dict[chrom].values())[0].prim_seed.cand_read_obj
        sread_name = sample_read.RNAME
        sread_seq = sample_read.seq
        stats_query = os.path.join(outdir, f"pis_params_sampleseq.que")
        with open(stats_query, "w") as que:
            que.write(f">{sread_name}_spanned_context\n")
            que.write(sread_seq + "\n")
        stats_out = os.path.join(outdir, f"pis_params_sampleseq.json")
        bl_params = get_blast_stats(stats_query, stats_out, genome, args)
    except:
        bl_params = "not available"
        ilog.vlprint(
            "INFO: Could not retrieve blast parameters in"
            " gen_integration_sites. This is not criticial.",
            2
            )

    pot_int_site_objs = list()
    for chrom in loc_dict:
        for nucl_site in loc_dict[chrom]:
            # use only the ReadToRefAl instances (tup[1]) instead of the whole  
            # tuple, since the coordinate and insertion-gap inormation is  
            # accesible in these objects via the self.int_coord and 
            # self.ins_gap attributes
            supp_als = loc_dict[chrom][nucl_site].all_als
            new_ins = PotIntSite(
                                chrom,
                                (
                                loc_dict[chrom][nucl_site].seed_coord_a, 
                                loc_dict[chrom][nucl_site].seed_coord_b
                                ),
                                supp_als,
                                bl_params,
                                args
                                )
            pot_int_site_objs.append(new_ins)
    
    ilog.vlprint(
        f"Generated {len(pot_int_site_objs)} potential integration site objects"
        , 3
    )

    return pot_int_site_objs


def merge_pot_int_sites(pot_int_site_objs, n=100, iteration=1):
    """Merge pot_int_sites objects that are less than n bp apart"""
    ilog.vlprint(f"Merge pot_int_sites objects mapping to same location. {iteration}.Iteration", 2)
    for site in pot_int_site_objs:
        for other in pot_int_site_objs:
            if site.chrom == other.chrom:
                if other != site:
                    try:
                        if site.approx_loc - n <= other.approx_loc <= site.approx_loc + n:
                            # check for evidence of multiple independent integrations in reads
                            #  belonging to PotIntSites poised for merging
                            multsite = any(
                                [al.cand_read_obj for al in site.supp_als if al.cand_read_obj.multiple_sep_int]
                                )
                            multother = any(
                                [al.cand_read_obj for al in other.supp_als if al.cand_read_obj.multiple_sep_int]
                                )
                            if any([multsite, multother]):
                                ilog.vlprint(
                                f"Omitted merging of potential integration sites {site.name} and "
                                f"{other.name}, because there was evidence for multiple independent"
                                f" integrations in reads associated with those sites", 4
                                )
                            else:
                                site.merge_with.add(other)
                                other.merge_with.add(site)
                    except TypeError:
                        ilog.vlprint(
                            "PotIntSites with less than 3 supporting reads cannot be merged.",
                            4, logging=False
                            )
    
    mergers = 0
    for site in pot_int_site_objs:
        if len(site.merge_with) > 0:
            site.merge(site.merge_with)
            mergers += 1
            site.merge_with = set()
    
    pot_int_site_objs = [site for site in pot_int_site_objs if not site.absorbed]
    ilog.vlprint(f"{mergers} site mergers were performed", 2)
    performed_mergers = True if mergers > 1 else False
    
    return pot_int_site_objs, performed_mergers


def filter_potintsites(pot_int_site_objs, min_supp, reads_info, args):
    """Filter all potential integration sites with regard to support quality"""
    ilog.vlprint("Filtering pot_int_site_objs", 2)
    before = len(pot_int_site_objs)
    supp_filt, loc_filt, imr_filt, avqual_filt  = 0, 0, 0, 0
    names = [pis.name for pis in pot_int_site_objs]
    pot_int_site_objs = [
        site for site in pot_int_site_objs if len(site.total_support) >= min_supp
        ]
    supp_filt += before - len(pot_int_site_objs)
    names2 = [pis.name for pis in pot_int_site_objs]
    ilog.vlprint(f"sites filtered out supp1: {[name for name in names if name not in names2]}", 6)
    for pint in pot_int_site_objs:
        pint.av_supp_qual = pint.eval_supp_qual()
    
    if not any([args.polyclonal, args.quick_mm2, args.intra]):
        piso_loc_filt = [pint for pint in pot_int_site_objs if pint.approx_loc]
        loc_filt += len(pot_int_site_objs) - len(piso_loc_filt)
        piso_imr_filt = [pint for pint in piso_loc_filt if pint.ident_match_ratio]
        imr_filt += len(piso_loc_filt) - len(piso_imr_filt)
        piso_avq_filt = [
            pint for pint in piso_imr_filt if pint.av_supp_qual > float(args.map_qual)
            ]
        avqual_filt += len(piso_imr_filt) - len(piso_avq_filt)
        pot_int_site_objs = [
            int_site for int_site in piso_avq_filt if len(int_site.supp_als) >= min_supp
            ]
        supp_filt += len(piso_avq_filt) - len(pot_int_site_objs)
    names3 = [pis.name for pis in pot_int_site_objs]
    ilog.vlprint(f"sites filtered out supp2: {[name for name in names if name not in names3]}", 6)

    dumped_reads = []
    for int_site in pot_int_site_objs:
        drs = int_site.dump_reads()
        dumped_reads.extend(drs)
    
    def parse_cand_csv(csv_file, header_start, col, sep=","):
        cands = []
        with open(csv_file, "r") as crv:
            for line in crv:
                if line.startswith(header_start):
                    pass
                else:
                    ls = line.strip().split(sep)
                    if len(ls) > 1:
                        if len("".join(ls)) > 1:
                            rname = ls[col]
                            cands.append(">"+rname+"\n")
        return cands
    
    if not args.quick_mm2 and not args.intra:
        all_cands = set(parse_cand_csv("candidate_reads", "RNAME", 2))
        filt_cands = set(parse_cand_csv("candidate_reads.fbr", "RNAME", 2))
        cont_genseq = set(parse_cand_csv("Integration_sites", "RNAME", 1))
    else:
        all_cands = set(parse_cand_csv("candidate_reads.paf", "", 0, sep="\t"))
        filt_cands = set(parse_cand_csv("candidate_reads_filt.tmp", "", 0, sep="\t"))
        cont_genseq = set(parse_cand_csv("Integration_sites.paf", "", 0, sep="\t"))


    not_dumped = [rname for rname in all_cands if rname not in dumped_reads]
    filt_for_len = [rname for rname in all_cands if rname not in filt_cands]
    nosup_passint = [rname for rname in not_dumped if rname not in filt_for_len]
    no_genseq = [rname for rname in nosup_passint if rname not in cont_genseq]
    nosuppassint2 = [rname for rname in nosup_passint if rname not in no_genseq]
    args.filt_for_len = [rn[1:].strip() for rn in filt_for_len]
    args.nosup_passint = [rn[1:].strip() for rn in nosuppassint2]
    args.no_genseq = [rn[1:].strip() for rn in no_genseq]

    ilog.vlprint(
        f"{len(dumped_reads)}/{len(all_cands)} candidate reads were written to" 
        f" site-specific supporting-reads files.", 4
        )
    reads_info["intread_supp_passedint"] = len(dumped_reads)
    ilog.vlprint(
        f"{len(filt_for_len)}/{len(all_cands)} candidate reads were filtered" 
        f" out because the length of int-sequences they contained were below"
        f" threshold.", 4
        )
    ilog.vlprint(
        f"{len(nosuppassint2)}/{len(all_cands)} candidate reads were written to" 
        f" 'candidate_reads_not_supporting_passed_int-site.fasta', because they do not"
        f" support integration sites that passed filters and evaluation. This"
        f" can have multiple reasons, including misalignment as well as"
        f" insufficient coverage for confirmation of an int-site.", 4
        )
    reads_info["intread_not_supp_passedint"] = len(nosuppassint2)
    ilog.vlprint(
        f"{len(no_genseq)}/{len(all_cands)} candidate reads were written to" 
        f" 'candidate_reads_without_sufficient_reference_seq.fasta', because"
        f" they cannot be mapped to the reference. This is usually due to"
        f" an insufficient portion of genomic sequences (e.g.reads consisting"
        f" exclusively of integration sequence).", 4
        )
    reads_info["intread_no_ref_map"] = len(no_genseq)

    if not args.quick_mm2 and not args.intra:
        get_reads_by_header(
            not_dumped,
            "Integration_candidate_reads.fasta",
            "Int-reads_not_dumped.fasta"
            )
        get_reads_by_header(
            filt_for_len,
            "Unfiltered_candidate_reads.fasta",
            "candidate_reads_filtered_for_short-int-cover.fasta"
            )
        get_reads_by_header(
            nosuppassint2,
            "Int-reads_not_dumped.fasta",
            "candidate_reads_not_supporting_passed_int-site.fasta"
            )
        get_reads_by_header(
            nosup_passint,
            "Int-reads_not_dumped.fasta",
            "candidate_reads_without_sufficient_reference_seq.fasta"
            )

    after = len(pot_int_site_objs)
    ilog.vlprint(f"{before - after} potential int-site/s were filtered out", 2)
    ilog.vlprint(
        f"Filter-type: supp={supp_filt}, loc={loc_filt}, imr={imr_filt},"
        f" avqual={avqual_filt}", 3
        )

    return pot_int_site_objs


def cleanup_aux_files(outdir, override, args):
    """Remove auxiliary files created during run not required anymore."""
    ilog.vlprint("Cleaning up auxiliary files.", 4)

    sub_dirs = ["logs", "figures", "reports", "reads_by_category"]
    sub_dirs = {sdir:"" for sdir in sub_dirs}
    for sdir in sub_dirs:
        sd_path = os.path.join(outdir, sdir)
        sub_dirs[sdir] = sd_path
        if not override:
            os.makedirs(sd_path, mode=777, exist_ok=False)
            os.chmod(sd_path, 0o777)

    aux_ext = {
        "log": sub_dirs["logs"],
        "html": sub_dirs["reports"],
        "pdf": sub_dirs["reports"],
        "csv": sub_dirs["reports"],
        "txt": sub_dirs["reports"],
        "json": sub_dirs["reports"],
        "svg": sub_dirs["figures"],
        "png": sub_dirs["figures"],
        "fasta": sub_dirs["reads_by_category"],
        }
    remove = [
        "tmp",
        "paf",
        "ava_res.csv",
        "Integration_sites",
        "prec_ints_pickled",
        "constr_chk.tmp",
        "constr_chk_refseq.out",
        "Cluster_not_assignable_reads.fasta",
        "ints_out.tmp",
        "IntSitesAlignments.json",
        "sec_site_blastdb_err.log",
        "Integration_sites_minimap2",
        "candidate_reads",
        "stats",
        "eval",
        "secsitedb",
        "cand_reads_without_integration",
        "ideo2png",
        "bait_candidate_reads",
        "pis_params_sampleseq",
        "Integration_candidate_reads.fasta",
        "Unfiltered_candidate_reads.fasta",
        "Int-reads_not_dumped.fasta",
    ]
    rem_empty = [
        "Multimer_info.txt",
        "multimeric_int.fasta",
        "multiple_indep_ints.fasta",
        "candidate_reads_not_supporting_passed_int-site.fasta",
        "candidate_reads_without_sufficient_reference_seq.fasta",
        "candidate_reads_filtered_for_short-int-cover.fasta"
    ]
    repl_except = {
        "Integration_summary_report.html": "",
        "single_chrom": "figures",
        "preex_single_chrom": "figures",
        "blast_map_figs": "figures",
        "minimap2_map_figs": "figures",
        "multi_intread_fig.html": "figures",
    }
    
    if not args.keep_nonint:
        remove.append("non_int_reads.fasta")

    for entry in os.scandir(outdir):
        fpath, tail = str(entry.path), os.path.split(str(entry.path))[1]
        tail_stem_start = [
            tail, tail.split(".")[0], tail.split(".")[-1], tail.split("_")[0]
            ]
        if tail in repl_except:
            os.replace(fpath, os.path.join(repl_except[tail], tail))
        elif any([word for word in tail_stem_start if word in remove]):
            os.remove(entry)
        elif tail in rem_empty and os.stat(fpath).st_size == 0:
            os.remove(entry)
        elif tail.split(".")[-1] in aux_ext:
            os.replace(fpath, os.path.join(aux_ext[tail.split(".")[-1]], tail))
        elif tail.split(".")[-1] == "tmp":
            os.remove(entry)
        else:
            pass


def check_retrieval_res(found_count):
    """Check if all searched reads could be retrieved successfully"""
    duplicate_count = {k:v for k,v in found_count.items() if v > 1}
    if len(duplicate_count) > 0:
        ilog.vlprint(
            "WARNING: Read file/s contain/s duplicate reads!!! "
            "Multiple use of identical read names may negatively affect the "
            "function of integration locator and generate (partially) incorrect"
            " results. Take care that every read has an individual, unique "
            "identifier, particularly when using reads from multiple sequencing"
            " runs, simulated reads or concatenated read files. The unique part"
            " of the read identifier must lie at the beginning of the header, "
            "before the first whitespace character."
            , 0
            )
    for entry in duplicate_count:
        ct = duplicate_count[entry]
        ilog.vlprint(f"WARNING: {ct} reads named {entry} were found!", 0)
    
    not_found = {k:v for k,v in found_count.items() if v == 0}
    if len(not_found) > 0:
        ilog.vlprint(
            "WARNING: Some reads could not be retrieved from read/s file/s."
            "If only a very minor fraction of candidate reads could not be "
            "retrieved this likely won't affect intloc results. If a substantial"
            "fraction of reads is missing better investigate."
            , 0
            )


def report_reads_info(
    reads_info, found_count_all, gen_ca_re_in, mmer_dct, full_int_len, trunc_int_len,
    prov_constr_len, args, opt_info={}
    ):
    cand_read_lens = [cr.slen for cr in gen_ca_re_in]
    av_crl = int(sum(cand_read_lens)/len(cand_read_lens)) if len(cand_read_lens)>0 else 0
    reads_info["av_len_intread"] = av_crl
    reads_info["int_read_filt_intlen"] = len(gen_ca_re_in)
    ilog.vlprint(
            f"A total of {len(found_count_all)} reads containing integration"
            f" sequences were detected.", 1
            )
    reads_info["int_read_all"] = len(found_count_all)
    ilog.vlprint(
            f"A total of {len(gen_ca_re_in)} integration site candidate reads"
            f" remain after filtering.", 1
            )
    ilog.vlprint(f"{mmer_dct['simple']} reads contain a simple integration", 1)
    ilog.vlprint(
        f"{len(mmer_dct['mult_ind'])} reads contain multiple separate integrations", 1
        )
    ilog.vlprint(
        f"{len(mmer_dct['multimer'])} reads contain a multimeric integration", 1
        )
    ilog.vlprint(
            f"Average length of integration candidate reads: {av_crl}", 1
            )
    ilog.vlprint(
            f"Length of provided integration-construct: {prov_constr_len}", 1
            )

    if args.int_len:
        av_int_len = args.int_len
        ilog.vlprint(
            f"Expected integration-length was manually set to {av_int_len}", 1
            )
        reads_info["av_len_int"] = av_int_len
    elif len(full_int_len) > 0:
        av_int_len = round(sum(full_int_len)/len(full_int_len))
        ilog.vlprint(
            f"Average length of (non-truncated) integrated construct in reads:"
            f" {av_int_len}", 1
            )
        reads_info["av_len_int"] = av_int_len
    else:
        ilog.vlprint(
            "No full representation of the integrating-construct could be found"
            " in the sequencing reads, possibly because all construct sequences"
            " are located at the end of reads. This may be expected if you are"
            " using a target enrichment strategy (e.g. CAS9 or PCR)."
            " Alternatively, the construct may become truncated during insertio"
            "n or length of reads in the sequencing library is not sufficient"
            " to span entire integrations. IntLoc will still be able to"
            " accurately place and identify integrations. If desired, for data"
            " where truncation of cunstruct is expected (e.g. enrichment) and"
            " homogenous int_len can be set manually.", 0
            )
        if len(trunc_int_len) > 0:
            av_int_len = round(sum(trunc_int_len)/len(trunc_int_len))
            ilog.vlprint(
            f"Average length of (truncated) integrated construct in reads:"
            f" {av_int_len}", 1
            )
            reads_info["av_len_int"] = av_int_len
        else:
            ilog.vlprint(
            "No part of the integrating-construct could be found"
            " in the sequencing reads.", 0
            )
            reads_info["av_len_int"] = "n.a."
            reads_info["av_len_intread"] = "n.a."
            av_int_len = -1
    reads_info["int_read_all"] = len(found_count_all)
    reads_info["int_read_filt_intlen"] = len(gen_ca_re_in)
    return av_int_len, av_crl


def exec_il_core(
    reads,
    construct,
    genome,
    outdir,
    min_supp,
    args,
    thr,
    reads_dir=False,
    bait=False
    ):
    """Run all processes required for identification of integrations"""
    tstart = report_time('Intloc core functions')

    def blast_n_retrieve(reads_path, construct, outdir, reads_dir=False, bait=bait):
        if bait:
            ilog.vlprint("blast and retrieve bait candidate reads", 1)
            bla_ca = blast_cand_reads(bait, reads_path, args, reads_dir=reads_dir, bait=bait)
            ca_fi = blast_output_filter(bla_ca[1], reads_dir=reads_dir, req_int_cont=args.ric)
            bait_res, filt_ct, all_ct = retrieve_cand_reads(
                bla_ca[1], ca_fi, reads_path, prefix="Bait"
                )
            return bait_res, filt_ct, all_ct
        else:
            ilog.vlprint("blast and retrieve integration candidate reads", 1)
            bla_ca = blast_cand_reads(construct, reads_path, args, reads_dir=reads_dir)
            ca_fi = blast_output_filter(bla_ca[1], reads_dir=reads_dir, req_int_cont=args.ric)
            ca_re_path, filt_ct, all_ct = retrieve_cand_reads(
                bla_ca[1], ca_fi, reads_path
                )
            return ca_fi, filt_ct, all_ct
    
    check_constr_for_refseqs(genome, construct, args, bait=bait)

    reads_info = {}
    found_count_all = Counter()
    found_count_filt = Counter()
    if reads_dir:
        fasta_ext = ["fasta", "fa", "fna", "fsa", "fas", "seq"]
        fastq_ext = ["fastq", "fq"]
        if len(list(os.scandir(reads_dir))) == 0:
            ilog.vlprint("The supplied read directory does not contain any files", 1)
            ilog.vlprint("Aborting integration locator.", 1)
            sys.exit()
        for entry in os.scandir(reads_dir):
            head, tail = os.path.split(str(entry.path))
            if tail.split(".")[-1] in fasta_ext and entry.is_file():
                ca_fi, flt_ct, all_ct = blast_n_retrieve(str(entry.path), construct, outdir, reads_dir=True)
                found_count_filt.update(flt_ct)
                found_count_all.update(all_ct)
            elif tail.split(".")[-1] in fastq_ext and entry.is_file():
                fq_conv = fastq_to_fasta(entry)
                fq_dir = os.path.join(head,"fastq") 
                if os.path.exists(fq_dir):
                    pass
                else:
                    os.makedirs(fq_dir, mode=777, exist_ok=False)
                    os.chmod(fq_dir, 0o777)
                os.replace(str(entry.path), os.path.join(fq_dir, tail))
                ca_fi, flt_ct, all_ct = blast_n_retrieve(fq_conv, construct, outdir, reads_dir=True)
                found_count_filt.update(flt_ct)
                found_count_all.update(all_ct)
            elif tail.split(".")[-1] == "gz" and tail.split(".")[-2] in fastq_ext and entry.is_file():
                comp_entry = entry
                entry = ".".join(str(entry.path).split(".")[:-1])
                with open(entry, "w") as extracted:
                    gio = gzip.open(comp_entry)
                    for line in gio:
                        extracted.write(line.decode("utf-8"))
                fq_conv = fastq_to_fasta(entry)
                if entry != str(comp_entry.path):
                    os.remove(entry)
                ca_fi, flt_ct, all_ct = blast_n_retrieve(fq_conv, construct, outdir, reads_dir=True)
                found_count_filt.update(flt_ct)
                found_count_all.update(all_ct)
            else:
                pass
        cleanup_blastdb(reads_dir)
    else:
        # This does not support .gz files yet; however these can be handled 
        # when the reads_dir arg (-d) is used instead of the file arg (-f)
        ca_fi, flt_ct, all_ct = blast_n_retrieve(reads, construct, outdir)
        found_count_filt.update(flt_ct)
        found_count_all.update(all_ct)
    
    check_retrieval_res(found_count_filt)
    
    if bait:
        cr_raw = sum(found_count_all.values())
        cr_filt = sum(found_count_filt.values())
        if cr_filt>0:
            ilog.vlprint(f"re-evaluate bait candidate reads with construct", 1)
            ca_fi, flt_ct, all_ct = blast_n_retrieve(
                ca_fi, construct, outdir, reads_dir=False, bait=False
                )
            cleanup_blastdb(outdir)
        else:
            ilog.vlprint(
                f"Bait based search for candidate reads did not yield any resul" 
                f"ts after filtering. Count before filtering: {cr_raw}. Exiting"
                f" intloc.", 0
                )
            sys.exit()

    gen_ca_re_in = gen_cand_read_instances(ca_fi, args)
    write_reads_wo_integration(gen_ca_re_in)
    mmers, mmer_dct = contains_multimers(gen_ca_re_in)

    full_int_len, trunc_int_len, prov_constr_len = write_intread_summary_txt(
        gen_ca_re_in, args
        )
    av_int_len, av_crl = report_reads_info(
        reads_info, found_count_all, gen_ca_re_in, mmer_dct, full_int_len, trunc_int_len, 
        prov_constr_len, args
        )

    outdir, cand_ints, genome = blast_int_genome(genome, args)
    int_site_als = gen_read_to_ref_al(cand_ints, gen_ca_re_in, args)

    loc_dict = nucleate_int_sites(int_site_als[1], args)
    pot_int_objs = gen_integration_sites(loc_dict, outdir, genome , args)

    if not args.polyclonal:
        perf_merge = True
        iterations = 0
        while perf_merge and iterations < 10:
            iterations += 1
            pot_int_objs, perf_merge = merge_pot_int_sites(
                pot_int_objs, n=args.min_dist_ints, iteration=iterations
                )
        pot_int_objs = filter_potintsites(
                pot_int_objs, min_supp, reads_info, args
                )
    else:
        import il_polyclonal
        il_polyclonal.ilog.verbosity = int(args.verbosity)
        pot_int_objs = filter_potintsites(
                pot_int_objs, min_supp, reads_info, args
                )
        pot_int_objs = il_polyclonal.run_il_polyclonal(
            genome, outdir, pot_int_objs, gen_ca_re_in, args
            )
        
    # currently write_output_txt() must be run before il_evaluate to generate 
    # IntSitesAlignments.json, which is used by load_alignments() in il_evaluate
    # to create IntSite objects; However, this could be ommited, if the information 
    # was handed to il_evaluate straight in the form of the PotIntSite-objects
    # without taking the detour through a file; this would also avoid to generate
    # conflicting output in the output ".txt" and other output formats
    # intcount_outtxt = write_output_txt(pot_int_objs)
    # int_dict = write_int_csv(pot_int_objs)
    write_output_txt(pot_int_objs)
    if args.polyclonal:
        ilog.vlprint(
        f"Skipping evaluation module for polyclonal sample.", 0
        )
    if len(pot_int_objs) < 1:
        ilog.vlprint(
        f" No potential integration sites found. Skipping evaluation module.", 0
        )
    if len(pot_int_objs) > 0 and not args.skip_eval:
        try:
            pot_int_objs = il_evaluate.run_il_evaluate(
                                                        outdir,
                                                        genome,
                                                        pot_int_objs,
                                                        gen_ca_re_in,
                                                        int_site_als[1],
                                                        reads,
                                                        av_int_len,
                                                        thr,
                                                        args
                                                        )
        except Exception as e:
            ilog.vlprint(f"WARNING: Error in il_evaluate {e}", 2)
            ilog.vlprint(
                f"WARNING: il_evaluate module failed. This can have multiple"
                f" reasons. Intloc will proceed with analysis and can still"
                f" generate valid data. However, some types of false positives"
                f" have a higher chance to evade elimination.", 2
                )
            for pio in pot_int_objs:
                pio.minimap_supp = "fail"
                pio.ident_match_ratio = "fail"
                pio.ssbr = "fail"
                pio.total_cov = "fail"
    else:
        for pio in pot_int_objs:
            if not args.polyclonal:
                pio.minimap_supp = "n.a."
            pio.ident_match_ratio = "n.a."
            pio.ssbr = "n.a."
            pio.total_cov = "n.a."

    int_cnt = len(pot_int_objs)
    int_dict = write_int_csv(pot_int_objs)
    write_cand_read_csv(gen_ca_re_in, args)
    draw_generic_ideogram(int_dict, pot_int_objs, args)
    try:
        for pio in pot_int_objs:
            pio.get_spanned_region_blast()
            # the aligner specific keys refer to the reads used to populate
            # the figure in il_map_figs, rather than the reads used to calculate
            # the boundaries of the span-reg (which are blast in both cases)
            ilog.vlprint(
                f"{pio.name}: spanned regions: minimap="
                f"{pio.spanned_region['minimap2']} blast="
                f"{pio.spanned_region['blast']}", 3
                )
        il_map_figs.gen_map_figs(pot_int_objs, aligner="blast")
    except:
        ilog.vlprint("INFO: Generation of blast map-figures failed", 1)
    tend = report_time('Intloc core functions', end=True)

    return tstart, tend, int_cnt, mmers, genome, pot_int_objs, av_crl, reads_info


def optional_modules(
    outdir,
    ava_clust,
    int_objs,
    av_crl,
    args,
    species=False,
    min_supp=3,
    testing=False,
    genome=False
    ):
    """Run optional modules depending on chosen options and supplied args"""
    features = None
    clust_objs = list()
  
    if ava_clust:
        report_time("All-vs-all read clustering")
        clust_objs = il_ava_clustering.ava_cluster(min_supp, args)
        clust_objs = il_ava_clustering.associate_cluster_with_potintsite(
            clust_objs, int_objs
            )
        report_time("All-vs-all read clustering", end=True)

    if species:
        import il_spec_ideo as ideo
        ideo.ilog.verbosity = int(args.verbosity)
        args.species = get_full_species_name(args)
        try:
            best_match, perfect = chk_assembly_ver(args.species, genome, args)
        except:
            best_match, perfect = False, False
            ilog.vlprint(
                "WARNING: Autoidentification and matching of reference assembly"
                " failed. That means no species specific ideogram can be"
                " generated or genomic features reported. Verification of"
                " integration positions in provided reference will proceed"
                " normally. See intlog.log for details.", 0
                )

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

        if best_match:
            import il_features as feat
            feat.ilog.verbosity = int(args.verbosity)
            generate_spec_ideogram(line_spacing, args)
        else:
            ilog.vlprint(
                "Species specific ideogram and feature table could not be"
                " generated, because the provided reference genome is not"
                " compatible with the available species resources.", 0
                )
    
    if testing:
        if not species:
            ilog.vlprint(
                "Species argument must be specified to perform testing. "
                "Testing skipped.", 2
                )
        else:
            import il_test
            il_test.ilog.verbosity = int(args.verbosity)
            try:
                il_test.compare_sim_ints_to_results(
                    args.species, outdir, testing, av_crl
                    )
            except UnboundLocalError as ule:
                ilog.vlprint(f"{ule}", 1)
                ilog.vlprint(
                    f"WARNING: Test routine failed, skipping test. This is not "
                    f"critical. This most likely happened due to missing csv-"
                    f"file with information about correct integration positions"
                    f". These files are generated by il_aux/il_interator.py,"
                    f" after artificialy inserting integrations into a"
                    f" reference sequence in silico. Per default such test-data"
                    f" is ony provided for yeast and human in intloc, but feel "
                    f"free to generate your own test-data for read-simulation "
                    f"and testing of intloc. Alternatively you can generate "
                    f"such a file manually for a known integration, by adhering"
                    f" to the correct csv format.  ", 0
                    )
 
    return len(clust_objs), args.species


def check_args(args):
    """Checking completness of file path arguments"""
    operating_system = platform.system()
    if operating_system == "Linux" or operating_system == "Darwin":
        path_sep = "/"
    elif operating_system == "Windows":
        path_sep = "\\"
    else:
        ilog.vlprint(
            "Integration locator has not been tested on the platform you are "
            "using, but may still work.", 0
            )
        path_sep = "/"

    path_args = ['reads_file', 'reads_dir', 'construct', 'genome', 'outdir',]
    for (key, val) in args.__dict__.items():
        if key in path_args:
            if val:
                if len(val.split(path_sep)) > 1:
                    pass
                elif len(val.split(path_sep)) == 0:
                    ilog.vlprint(
                        f"Path to required argument {key} is missing. Try to "
                        f"supply full path to file.", 0
                        )
                    ilog.vlprint("Aborting program", 0)
                    sys.exit()
                else:
                    args.__dict__[key] = os.path.join(args.run_dir, val)


def run_intloc(args):
    """integration_locator main function"""
    
    pkg_dir = os.path.split(os.path.abspath(__file__))[0]
    check_args(args)
    outdir = os.path.join(
        args.outdir, "intloc_out_" + strftime("%Y-%m-%d_%H-%M-%S", localtime())
        )
    args.outdir = outdir
    verbosity = int(args.verbosity)
    res_dir = os.path.join(pkg_dir, "il_resources")
    args.res_dir = res_dir
    args.pkg_dir = pkg_dir
    args.defaults = DEFAULTS
    initialize_loggers(args) # does not initialize the logger of the main module
    override, gui, bait, = args.override, args.gui, args.bait
    cleanup = args.cleanup_aux
    # could the next line be omitted if integration_locator module was added to 
    # log_init_modules.json?
    ilog.outdir, ilog.verbosity = outdir, verbosity
    check = args.check_reads
    thr = determine_num_threads(args)

    if args.reads_file:
        reads = args.reads_file
        gen_output_dir(outdir, reads_dir=False, override=override, gui=gui)
        os.chdir(outdir)
        reads = test_input(reads, check)
    else:
        reads = args.reads_dir
        gen_output_dir(outdir, reads_dir=True, override=override, gui=gui)
    os.chdir(outdir)
    ilog.vlprint(f"intloc version {__version__}", 1)
    ilog.vlprint(f"intloc command line: {' '.join(sys.argv)}", 10)
    req_deps, fac_deps = required_dependencies(args)
    for rdep in req_deps:
        chk_nonpy_dependency(rdep, args)
    for fdep in fac_deps:
        chk_nonpy_dependency(fdep, args)

    bait, construct = args.bait, args.construct
    args.construct = test_input(construct, check, constr_fa=True)
    construct = args.construct
    min_supp, Nx, = args.min_supp, float(args.Nx)
    if args.polyclonal:
        min_supp = 1
        args.skip_eval = True
        args.min_dist_ints = 10
        args.sample = "polyclonal"
    else:
        args.sample = "clonal"
    genome, dwnl_genome = args.genome, args.dwnl_genome,
    if dwnl_genome:
        dwnl_success, genome = run_getasm(dwnl_genome)
        args.genome = genome
    if args.coverage:
        if args.coverage > 15:
            min_supp = int(args.coverage / 5)
    species, ava_clust = args.species, args.cluster_reads
    cleanup, testing = args.cleanup_aux, args.testing

    if args.quick_mm2:
        chk_compatibility_args(
            ["intra", "ava_clust", "polyclonal", "reads_dir"], "quick", args,
            opt_msg=[
                "Reads must be supplied as single fasta file.",
                "If you want to analyse polyclonal data in fast mode, set --min"\
                "_read_supp to 1 and --min_dist_ints to a reasonably low value"\
                " <50 (e.g. 10)."
                ]
            )
        set_module_defaults("quick", ["min_dist_ints"], args)
        import il_quick_mm2 as iq
        iq.ilog.verbosity = int(args.verbosity)
        qui_res = iq.run_quick_mm2(
                                reads,
                                construct,
                                genome,
                                outdir,
                                min_supp,
                                args,
                                )
        int_objs, crl, start, end, intcnt, multmers, readinf = qui_res
    elif args.intra:
        import il_intra
        il_intra.ilog.verbosity = int(args.verbosity)
        chk_compatibility_args(
            ["quick", "ava_clust", "polyclonal", "reads_dir"], "intra",
            args, opt_msg=["Reads must be supplied as single fasta file."]
            )
        args.skip_eval, args.color_ints = True, "lime"
        set_module_defaults("intra", ["min_dist_ints", "req_int_cont"], args)
        int_res = il_intra.run_il_intra(args)
        int_objs, crl, start, end, intcnt, multmers, readinf, opt_info = int_res
    else:
        cor_res = exec_il_core(
                            reads,
                            construct,
                            genome,
                            outdir,
                            min_supp,
                            args,
                            thr,
                            reads_dir=args.reads_dir,
                            bait=bait
                            )
        start, end, intcnt, multmers, genome, int_objs, crl, readinf = cor_res
    
    clust, species = optional_modules(
                                    outdir,
                                    ava_clust,
                                    int_objs,
                                    crl,
                                    args,
                                    species=species,
                                    min_supp=min_supp,
                                    testing=testing,
                                    genome=genome
                                    )

    if not args.intra:
        opt_info = {}

    write_run_info(
            start,
            end,
            reads,
            construct,
            genome,
            outdir,
            min_supp,
            intcnt,
            clust,
            multmers,
            readinf,
            thr,
            args,
            opt_info=opt_info
            )

    try:
        write_pdf(Nx, species=species, avaclust=ava_clust)
    except Exception as e:
        ilog.vlprint("WARNING: PDF summary report could not be created", 0)
        ilog.vlprint(str(e), 0)

    try:
        il_html.write_html_report(args)
    except Exception as e:
        ilog.vlprint("WARNING: html summary report could not be created", 0)
        ilog.vlprint(str(e), 0)
    
    if args.gene_dens:
        try:
            il_mod_ideo.run_mod_ideo(args)
        except:
            ilog.vlprint("Derivation of gene-density figures failed.", 1)
    
    if args.feat_dens:
        try:
            targetfeat, targetclass = args.feat_dens.split(",")
            il_mod_ideo.run_mod_ideo(args, tfeat=targetfeat,tclass=targetclass)
        except:
            ilog.vlprint("Derivation of feat-density figures failed.", 1)

    defs = args.__dict__.pop("defaults")
    ilog.vlprint(f"intloc defaults: {defs}", 10)
    ilog.vlprint(f"intloc run parameters: {args}", 10)
    report_time("preparation of summary reports")
    if cleanup:
        cleanup_aux_files(outdir, override, args)
        log = os.path.join("logs", "intlog.log")
        report_time("Preparation of summary reports",end=True, alt_log=log)
        ilog.vlprint(f"Results saved to: {args.outdir}", 0, alt_log=log)
    else:
        report_time("Preparation of summary reports", end=True)
        ilog.vlprint(f"Results saved to: {args.outdir}", 0)
    
    if args.cat_log:
        with open(os.path.join(args.outdir, "logs", "intlog.log"), "r") as log:
            data = log.read()
            print(data)


def main():
    "Integration locator main function"
    args = get_arguments()
    run_intloc(args)


if __name__ == "__main__":
    main()
    
