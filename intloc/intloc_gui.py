#!/usr/bin/env python3

# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import argparse
import sys
import os
from gooey import Gooey, GooeyParser

# import il_evaluate
sys.path.append(os.path.split(__file__)[0])
from version import __version__
from integration_locator import *
from il_chk_conv import ranged_typechk


img_dir = os.path.join(os.path.split(__file__)[0], 'il_resources')

@Gooey(
    advanced=True,
    program_name='Integration Locator',
    # navigation='TABBED', # only for subparsers
    tabbed_groups=True, # for argument groups
    image_dir= img_dir,
    richtext_controls=True,
    menu=[ 
        {
        "name": "IntLoc About",
        "items": [
                {
                "name": "IntLoc",
                "type": "AboutDialog",
                "menuTitle": "IntLoc",
                "description": "IntLoc the definitive integration locator",
                "version": __version__,
                "copyright": '2022',
                "website": 'https://github.com/simakro/Integration_locator',
                "developer": "Simon Magin",
                "license": 'BSD2'
                },
                {
                "name": "Help",
                "menuTitle": "Help",
                "type": "Link",
                "url": 'https://github.com/simakro/Integration_locator',
                }
            ]
        }
    ]
    )
def get_arguments():
    parser = GooeyParser(
        description="Integration locator detects the number and location of "
        "integrations of transgenic constructs or other inserted DNA sequences "
        "within a genome using long-read sequencing data.",
        add_help=True
        )

    required_args = parser.add_argument_group('Required_arguments')
    figure_args = parser.add_argument_group('Figures')
    settings = parser.add_argument_group('Settings')
    optional_args = parser.add_argument_group('Optional_arguments')
    modes = parser.add_argument_group('IntLoc_Modes')

    genome_choice = required_args.add_mutually_exclusive_group(required=True)
    genome_choice.add_argument(
        "-g", "--genome", help="Path to reference genome",
        widget="FileChooser"
        )
    genome_choice.add_argument(
        "-w", "--dwnl_genome",
        help="Specify a species for download of latest assembly", 
        widget="Dropdown",
        choices=[
            "Homo sapiens",
            "Mus musculus",
            "Saccharomyces cerevisiae",
            "Arabidopsis thaliana",
            "Rattus norvegicus",
            "Drosophila melanogaster",
            "Danio rerio",
            "Caenorhabditis elegans"
            ]
         )
    
    read_input = required_args.add_mutually_exclusive_group(required=True)
    read_input.add_argument(
        "-f", "--reads_file", help="Path to a single fasta/q file containing"
        " all sequencing reads", widget="FileChooser"
        )
    read_input.add_argument(
        "-d", "--reads_dir", help="Path to a directory containing multiple"
        " fasta/q files with reads. Not recursive, subdirectories are ignored."
        " (mutually exclusive with -f/--reads_file)", widget="DirChooser"
         )

    required_args.add_argument(
        "-c", "--construct", required=True,
        help="Path to fasta file containing the sequence  of the integrating "
        "construct", widget="FileChooser"
        )
    required_args.add_argument(
        "-o", "--outdir", required=True, help="Path for output directory",
        widget="DirChooser"
        )
    # min_read_supp arg is optional in CLI version, but required in intloc_gui 
    required_args.add_argument(
        "--min_read_supp", dest='min_supp', default=DEFAULTS["min_supp"],
        help=f"Set minimum read support for integration validation"
        f" [{DEFAULTS['min_supp']}]", type=ranged_typechk(1,1000,int)
        )

    figure_args.add_argument(
        "--spec", dest='species', default=None, widget="Dropdown",
        choices=[
            "human",
            "mouse",
            "yeast",
            "arabidopsis",
            "rat",
            "drosophila",
            "zebrafish",
            "c.elegans"
            ],
        help="Select species (e.g. \"human\") to enable output of species "
        "specific ideogram and annotations"
        )
    figure_args.add_argument(
        "--color_ints", dest='color_ints', default="lime",widget="Dropdown",
        choices=[
            "lime",
            "green",
            "blue",
            "cyan",
            "magenta",
            "yellow",
            "olive",
            "orange",
            "pink",
            "purple",
            "teal",
            "violet"
            ],
        help="Specify color of int labels in figures. Not available in intra-mo"
        "de.[default='lime']"
        )
    figure_args.add_argument(
        "--gene_dens", dest='gene_dens', default=False, action='store_true',
        help="Derive ideograms with gene-density information."
        )
    figure_args.add_argument(
        "--feat_dens", dest='feat_dens', default=False,
        type=boolmix_typechk(str),
        help="Derive ideograms with feature-density information. A feature clas"
        "s and subclass have to be provided separated by a comma. (--feat_dens "
        "gene,protein_coding reproduces the function of the --gene_dens flag)"
        )
    figure_args.add_argument(
        "--ntigs", dest='ntigs', default=False, 
        type=ranged_typechk(0,1000,int, alt=bool),
        help="Set number of contigs to be included in generic ideogram. If not"
        " set this number is handled by the Nx parameter. [Default=False]"
        )
    figure_args.add_argument(
        "--include_Nx", dest='Nx', default=98,
        widget="Slider", type=ranged_typechk(1,100,float),
        help="Choose cutoff for contigs to be included in figures in generic id"
        "eogram. Contigs comprising at least x percent of bases in RefGenome "
        "[98] will be included. This parameter can be helpful to avoid graphs "
        "becoming overcrowded when working with reference genomes that contain "
        "large amounts of alternative loci or a high number of miniscule "
        "unplaced contigs. Don't worry, this applies only to the figures 1-3 in"
        " the pdf-report. Thus, no sequences will be ignored during integration"
        " search, i.e. there is no risk of missing out on integration sites. "
        "Species specific ideogram is not affected."
        )

    settings.add_argument(
        "--verbosity", dest='verbosity', default=1,
        type=ranged_typechk(-1,100,int), 
        help="Adjust verbosity of output during runtime"
        )
    settings.add_argument(
        "--file_check", dest='check_reads', default="quick",
        help="Adjust thoroughness of file-check; quick/extensive [quick]"
        )
    settings.add_argument(
        "--skip_eval", dest='skip_eval', default=False, action='store_true',
        help="Skip evaluation module. Evluation for the identification of false"
        " positives can significantly increase the runtime (in dependence of"
        " genome size and number of integrations). You can skip it when in a"
        " hurry and do not fear no FPs.", gooey_options={'visible': False}
        )
    settings.add_argument(
        "--cores", dest='cores', default=False,
        type=ranged_typechk(1, 1024, int, out_type=str, alt=bool),
        help="Specify number of cores to use. Intloc by default attempts to"
        " run parallelized tasks on all but two of the cores available on the"   
        " current machine (e.g. using 6 of 8 available cpu cores)."
        )
    settings.add_argument(
        "--cleanup", dest='cleanup_aux', default=True, action='store_true',
        help="Cleanup auxiliary files [True]. Relevant for debugging only.", 
        gooey_options={'visible': False}
        )


    optional_args.add_argument(
        "--bait", dest='bait', widget="FileChooser",
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
        "-s", "--seq_tec", dest='seq_tec', default="ont", widget="Dropdown",
        choices=["ont", "pb"],
        help="Specify (long-read) sequencing technology (ont/pb). "
        "[default='ont']"
        )
    optional_args.add_argument(
        "--min_dist_ints", dest='min_dist_ints', default=25,
        type=ranged_typechk(25,50000,int),
        help="Sets the minimimun distance in bp that is required in order to "
        " count integrations on the same DNA molecule as separate entities or"
        " part of a complex integration event. [defaults by mode: standard=25,"
        " quick=1000]"
        )
    optional_args.add_argument(
        "--int_len", dest='int_len', default=False,
        type=ranged_typechk(25,10000000,int, alt=bool),
        help="Manually set the expected length of integrations. For WGS data"
        " this is not required nor recommended, since IntLoc automatically"
        " determines the lengths of integrations in the genome. However, this"
        " parameter can be helpful when only using part of the integration as"
        " bait, or when using enrichment methods (e.g. Cas9) that always result"
        " in truncated integration sequences in the reads. [default=False]"
        )
    optional_args.add_argument(
        "--untrimmed", dest='untrimmed', default=0,
        type=ranged_typechk(0,250,int),
        help="If you are using untrimmed raw-reads as input, you can provide "
        " the expected lengths of sequencing adapters + barcodes (if used) at"
        " the ends of each read via this parameter. [default=0]"
        )
    optional_args.add_argument(
        "-q", "--map_qual", dest='map_qual', default=90,
        help="Specify cutoff for required completeness (not identity) of"
        " mapping. If the genome of your sample does not deviate significantly"
        " from the reference, this value should be 3-10 %% lower than the"
        " accuracy of the sequencing reads. Increasing the tolerance by" 
        " lowering this value may help to reveal integrations in genomic"
        " locations of your sample that significantly deviate from the"
        " reference. However, it may increase the risk for false positives as"
        " well. [default=90]"
        )
    optional_args.add_argument(
        "--no_multi_ints_info", dest='no_multi_ints', default=False,
        action='store_true',
        help="Toggle precise analysis of reads with multiple insertions of the "
        "integrating sequence. When deactivated those reads will typically just"
        " support one integration site and linkage analysis does not take place"
        ". [default=True]"
        )
    optional_args.add_argument(
        "--ava_clust", dest='cluster_reads', default=False, action='store_true', 
        help="Cluster reads based on all-vs-all alignment. Comparison of read-"
        "clusters to integration-sites identified by alignment of reads to the "
        "reference genome may assist identification of potential duplicates "
        "that could arise from integration into highly repetetive genomic loci."
        )
    optional_args.add_argument(
        "--req_int_cont", dest='ric', default=20,type=ranged_typechk(1,100,int),
        help="Only include sequencing reads containing at least x%% of the int"
        "egrating sequence as candidates for the detection of integration sites"
        ". [default=20]"
        )
    optional_args.add_argument(
        "--cov", dest='coverage',
        default=None, type=ranged_typechk(1,100000, float),
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

    il_mode = modes.add_mutually_exclusive_group()
    il_mode.add_argument(
        "--standard", dest='standard', default=False, action='store_true',
        help="IntLoc standard mode, for identification of transgenic integratio"
        "ns in monoclonal samples (WGS or target-enriched). [Default]"
        )
    il_mode.add_argument(
        "--polyclonal", dest='polyclonal', default=False, action='store_true',
        help="Polyclonal mode. Use to identify integrations in DNA sampled from"
        " a heterogenous mixture of cells. (Sets min_read_supp=1, prevents merg"
        "ing of int-sites that are spaced less than min_dist_ints from each oth"
        "er and prevents exclusion of sites that share extensive similarity and"
        " overlap."
        )
    il_mode.add_argument(
        "--intra", dest='intra', default=False, action='store_true',
        help="Intra-mode allows to find pre-existing as well as newly occuring"
        " integrations of a provided query sequence in the genome. This mode ca"
        "n for example be used to identify the locations of LINEs in a genome."
        )
    il_mode.add_argument(
        "--quick", dest='quick_mm2', default=False, action='store_true',
        help="Quick and dirty integration detection (using only minimap2) witho"
        "ut thorough evaluation of sites. Often sufficient for easy cases (smal"
        "ler non-repetetive genomes, low number of integrations), but typically"
        " produces poorer results in more demanding cases (large complex genome"
        "s, high number of integrations) than standard mode."
        "Not compatible with following arguments: --read_dir, --ava_clust"
        )
    
    optional_args.add_argument(
        "--override", dest='override', default=False, action='store_false', 
        help=argparse.SUPPRESS, gooey_options={'visible': False}
        )
    optional_args.add_argument(
        "--no_culling", dest='no_culling', default=False, action='store_true',
        help=argparse.SUPPRESS,  gooey_options={'visible': False}
        )
    optional_args.add_argument(
        "--unfiltered", dest='unfiltered', default=False, action='store_true',
        help=argparse.SUPPRESS, gooey_options={'visible': False}
        )
    optional_args.add_argument(
        "--gui", dest='gui', default=True, action='store_true', 
        help=argparse.SUPPRESS, gooey_options={'visible': False}
        )
    optional_args.add_argument(
        "--testing", dest='testing', default=False, action='store_true', 
        help=argparse.SUPPRESS, gooey_options={'visible': False}
        )
    optional_args.add_argument(
        "--cat_log", dest='cat_log', default=False, action="store_true",
        # help="Print intlog.log to stdout"
        help=argparse.SUPPRESS, gooey_options={'visible': False}
        )
    optional_args.add_argument(
        "--dep_ver", dest='dep_ver', default={},
        # help="dict for storage of infos about dependency versions"
        help=argparse.SUPPRESS, gooey_options={'visible': False}
        )

    args = parser.parse_args()

    return args


def gui_main():
    "Integration locator main function"
    args = get_arguments()
    run_intloc(args)


if __name__ == "__main__":
    #pkg_dir = os.path.split(os.path.abspath(__file__))[0]
    # args = get_arguments()
    # run_intloc(args)
    gui_main()
