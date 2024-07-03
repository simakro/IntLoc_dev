# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import platform
import os
import math
import subprocess
import matplotlib.pyplot as plt

from il_logging import Ilogger
from version import __version__
# try:
    # from intloc.version import __version__
# except ModuleNotFoundError:


ilog = Ilogger()
ilog.module = __name__

def write_int_csv(pot_int_site_objs, prefix=""):
    """Generate csv integration report"""

    def eval_typechk(mms_val):
        if type(mms_val) == dict:
            return len(mms_val)
        elif type(mms_val) == str:
            return mms_val
        elif type(mms_val) == int:
            return str(mms_val)
        elif type(mms_val) == float:
            return str(mms_val)
        elif type(mms_val) == bool:
            if not mms_val:
                return "failed"
            else:
                return "unexpected value"
        else:
            return "type unexpected"

    ilog.vlprint("Writing integration_report.csv and creating int_dict", 3)
    with open(f'{prefix}Integration_Report.csv', 'w') as rep:
        rep.write(
            "ID,INT_name,CHR/contig,loc,+/-[bp],int-rs,rs(m),non-int-r,ts_cov,IMR,SSBR,read_names\n"
            )
        int_contigs = []
        id = 0
        for int_site in pot_int_site_objs:
            id += 1
            if len(int_site.supp_als) > 0:
                supp_reads = [al.READ for al in int_site.supp_als]
            else:
                supp_reads = [al for al in int_site.minimap_supp]
            rep.write(
                f"{id},"
                f"{int_site.name}," 
                f"{int_site.chrom}," 
                f"{int_site.approx_loc},"
                f"{int_site.app_loc_sd},"
                f"{len(int_site.supp_als)},"
                f"{eval_typechk(int_site.minimap_supp)},"
                f"{int_site.sum_non_int_als},"
                f"{int_site.total_cov},"
                f"{str(int_site.ident_match_ratio)[:5]},"
                f"{str(eval_typechk(int_site.ssbr))[:5]},"
                f"{(';'.join(supp_reads))}\n"
                )
            int_contigs.append(int_site.chrom)
        set_contigs = set(int_contigs)
        int_dict = {
            chrom: int_contigs.count(chrom) for chrom in sorted(set_contigs)
            }
    return int_dict

def write_cand_read_csv(cand_read_lst, args):
    """Generate csv report all candidate reads carrying integrated sequence"""
    ilog.vlprint(
        "Generating csv report for all candidate reads carrying integrated sequence", 3
        )
    with open('Candidate_Reads.csv', 'w') as rep:
        rep.write(
            "Read_ID, Read_name, read_length, int_pos_in_read, int_len_in_read,"
            " int_cov_bases, int_ident, evalue, bitscore, qcovhsp, qcovus,"
            " completeness\n"
            )
        for al in cand_read_lst:
            al_sum = al.report_al_sum()
            # print(al.RNAME, al.read_cov_range)
            int_pos_in_read = str(al_sum[3][0]) + "-" + str(al_sum[3][1])
            # print("int_pos_in_read", int_pos_in_read)
            int_cov_bases = str(al_sum[5][0]) + "-" + str(al_sum[5][1])
            completeness = al.trunc_constr_read_pos()
            rep.write(
                f"{al_sum[0]}, {al_sum[1]}, {al_sum[2]}, {int_pos_in_read},"
                f" {al_sum[4]}, {int_cov_bases}, {al_sum[6]}, {al_sum[9]},"
                f" {al_sum[10]}, {al_sum[11]}, {al_sum[12]}, {completeness}\n"
                )


def write_run_info(
    start_time, end_time, reads, construct, genome, outdir, min_supp, int_count,
    clust_count, multimers, reads_info, thr, args, opt_info={}
    ):
    """Write run info csv"""
    ilog.vlprint("Writing run info", 3)
    modes = {"intra": args.intra, "polyclonal": args.polyclonal, "quick": args.quick_mm2}
    modes = [k for k,v in  modes.items() if v]
    mode = "standard" if len(modes)==0 else " ".join(modes)
    try:
        blast_ver = subprocess.Popen([args.blastn, "-version"], stdout=subprocess.PIPE)
        get_bv = blast_ver.stdout.read().replace(b',', b';').replace(b'\n', b'').replace(b'\r', b'').decode("utf-8")
    except:
        get_bv = "n.a."
    try:
        mm2_ver = subprocess.Popen([args.minimap2, "--version"], stdout=subprocess.PIPE)
        get_mm = mm2_ver.stdout.read().replace(b',', b';').replace(b'\n', b'').replace(b'\r', b'').decode("utf-8")
    except:
        get_mm = "n.a."
    with open("run_info.csv", "w") as inf:
        inf.write(f"nodename,{platform.uname().node}\n")
        inf.write(f"OS,{platform.uname().system}\n")
        inf.write(f"sys_release,{platform.uname().release}\n")
        inf.write(f"sys_version,{platform.uname().version}\n")
        inf.write(f"IntLoc_version,{__version__}\n")
        inf.write(f"python_ver,{platform.python_version()}\n")
        inf.write(f"BLAST_ver,{get_bv}\n")
        inf.write(f"minimap2_ver,{get_mm}\n")
        inf.write(f"num_CPUs,{os.cpu_count()}\n")
        inf.write(f"machine,{platform.uname().machine}\n")
        inf.write(f"processor,{platform.uname().processor}\n")
        inf.write(f"cores_used,{thr}\n")
        inf.write(f"start_time,{start_time[0]}\n")
        inf.write(f"end_time,{end_time[0]}\n")
        inf.write(f"run_time (sec),{end_time[1] - start_time[1]}\n")
        inf.write(f"mode,{mode}\n")
        inf.write(f"Sequencing Reads,{reads}\n")
        inf.write(f"Integrating Sequence,{construct}\n")
        if args.bait:
            inf.write(f"Bait Sequence,{args.bait}\n")
        inf.write(f"Reference Genome,{genome}\n")
        inf.write(f"Output directory,{outdir}\n")
        inf.write(
            f"Number of reads containing evidence of integration sequence,"
            f"{reads_info['int_read_all']}\n"
            )
        if not args.intra:
            inf.write(
                f"Number of int-reads passing integration-sequence length filter,"
                f"{reads_info['int_read_filt_intlen']}\n"
                )
            inf.write(
                f"Int-reads not mapping to reference,"
                f"{reads_info['intread_no_ref_map']}\n"
                )
            inf.write(
                f"Other int-reads not supporting a confirmed Integration,"
                f"{reads_info['intread_not_supp_passedint']}\n"
                )
            inf.write(
                f"Int-reads supporting confirmed Integration sites,"
                f"{reads_info['intread_supp_passedint']}\n"
                ) 
        inf.write(
            f"Average length of (non-truncated) integrated construct in reads,"
            f"{reads_info['av_len_int']}\n"
            )
        inf.write(
            f"Average length of integration candidate reads,"
            f"{reads_info['av_len_intread']}\n"
            )
        inf.write(f"Minimum read support,{min_supp}\n")
        if args.intra:
            # inf.write(
            # f"Number of Integrations in reference, {opt_info['Integrations in genomic reference']}\n"
            # )
            inf.write(
            f"Number of novel Integrations,{int_count}\n"
            )
        else:
            inf.write(
                f"Number of confirmed Integrations,{int_count}\n"
                )
        if clust_count:
            inf.write(f"# Read Clusters (ava-alignm.),{clust_count}\n")
        if not any([args.intra, args.quick_mm2]):
            if multimers:
                inf.write(
                    "Multimer Info,Indication for multimeric or merged "
                    "integrations exists.\n"
                    )
            else:
                inf.write(
                    "Multimer Info,No evidence for multimeric or merged "
                    "integrations detected.\n"
                    )
        for key in opt_info:
            inf.write(
                    f"{key},{opt_info[key]}\n"
                    )


def get_genome_info_blast(genome_db):
    """Get reference info from BLAST DB"""
    ilog.vlprint("Retrieving genome info from BLAST DB", 3)
    # Collect infos about RefGenome database formatted to contain 
    # SeqID, SeqTitle, SeqLength, SeqHashval ("%i %t %l %h")
    db_entr = subprocess.Popen(
        [
            "blastdbcmd",
            "-db",genome_db,
            "-entry", "all",
            "-outfmt", "%i::$&~::%t::$&~::%l::$&~::%h"
        ],
        stdout=subprocess.PIPE
        )
    dbe_read = db_entr.stdout.readlines()
    entr_num = 0
    entries_dict = {}
    # extract Contig name and length and save in entries dict
    for entry in dbe_read:
        entry = entry.decode("utf-8")
        entr_num += 1
        espl = entry.split("::$&~::")
    # gnl = general db identifier; lcl= local seq identifier
    # ref= NCBI reference sequence; bbs = GenInfo Backbone Id
    # for GenBank, EMBL Data Library and DDBJ:
    # gi = gi|gi-number|gb or emb or dbj|accession|locus
        if entry.startswith("gnl"):
            cont_id, descr =  espl[1], espl[2]
        else:
            cont_id, descr = espl[0], espl[1]
        cont_len, cont_hash = espl[-2], espl[-1]
        # Convert blastdbcmd title output from newer versions (e.g. 2.11.0; 
        # "ref|name|) or older versions (e.g. <2.9;
        # "lcl|name| to just contig name
        if cont_id.startswith("ref"):
            cont_id = cont_id.split("|")[1]
        if cont_id.startswith("lcl"):
            cont_id = cont_id.split("|")[1]
        if cont_id.startswith("gb"):
            cont_id = cont_id.split("|")[1]
        if cont_id.startswith("ENA"):
            cont_id = cont_id.split("|")[2]
        # the below two have not yet been verified
        if cont_id.startswith("emb"):
            cont_id = cont_id.split("|")[1]
        if cont_id.startswith("dbj"):
            cont_id = cont_id.split("|")[1]
        entries_dict[cont_id] = {
                                "descr": descr,
                                "len": int(cont_len),
                                "hash": cont_hash,
                                "source": "blastdbcmd",
                                }
    return entries_dict


def scan_ref_genome(genome_db, int_dict):
    """Parse reference genome to retrieve contig names and lengths."""
    ilog.vlprint(
        "Parsing reference genome to retrieve contig names and lengths", 1
        )
    entries_dict = {}
    entries_hash = {}
    curr_cont = ""
    with open(genome_db, "r") as gen:
        for line in gen:
            if line.startswith(">"):
                curr_cont = line.strip()
                entries_dict[curr_cont] = 0
                entries_hash[curr_cont] = []
            else:
                entries_dict[curr_cont] += len(line.strip())
                entries_hash[curr_cont].append(line.strip())
    entries_hash = {k:"".join(v) for k,v in entries_hash.items()}
    entries_hash = {k:hash(v) for k,v in entries_hash.items()}

    idx_saccver = []
    for contig in entries_dict:
        contig_lst = contig.replace(">", "").split(" ")
        for _ in contig_lst:
            if _ in int_dict:
                idx_saccver.append(contig_lst.index(_))

    new_entries = {}
    count_idx = {idx_saccver.count(n): n for n in idx_saccver}
    if len(count_idx)>0:
        max_idx = count_idx[max(count_idx)]
    else:
        max_idx = 0
    for contig in entries_dict:
        split_tig = contig.replace(">", "").split(" ")
        new_name = str(split_tig[max_idx])
        split_tig.pop(max_idx)
        descr = " ".join(split_tig)
        new_entries[new_name] = {
                                "len": entries_dict[contig],
                                "descr": descr,
                                "hash": entries_hash[contig],
                                "source": "scan_ref_db()",
                                }

    stats = [entries_dict[tig] for tig in entries_dict]
    total_bp = sum(stats)
    ilog.vlprint("Genome reference statistics:", 10)
    ilog.vlprint(f"No. of contigs {len(stats)}", 10)
    ilog.vlprint(f"longest contig {max(stats)}", 10)
    ilog.vlprint(f"Total base count {total_bp}", 10)
    return new_entries


def select_contigs(entries_dict, int_dict, x=98, n=False):
    """Reduce reference information to a set of contigs for graphic  
    representation that comprise together at least x% [default 98%] of the 
    genome and ignore small unplaced fragments, unless they carry a potential 
    integration."""
    # Calculate total size (bp) of Genome database
    select_tigs = {}
    sum_db_size = 0
    for k in entries_dict:
        sum_db_size += entries_dict[k]["len"]
    sorted_entries = sorted(entries_dict.items(), key=lambda n: n[1]["len"])

    if n:
        select_entries = sorted_entries[-n:]
        select_tigs = {k:v["len"] for k,v in select_entries}
    else:
        ilog.vlprint(
            "Selecting set of reference contigs according to Nx parameter for "
            "graphical representation", 5
            )
        # Include contigs until the selected contigs comprise >= x % of the genome 
        # => Exclude potentially large number of tiny unplaced fragments
        # select_tigs = {}
        sorted_entries = sorted(entries_dict.items(), key=lambda n: n[1]["len"])
        acc_len = 0
        for i in range(len(sorted_entries)-1, 0, -1):
            if acc_len < (sum_db_size/100)*x:
                acc_len += sorted_entries[i][1]["len"]
                select_tigs[sorted_entries[i][0]] = sorted_entries[i][1]["len"]
            else:
                break
        # Should any of the excluded fragments carry an integration, include it
        #  into selection as well
    for int_contig in int_dict:
        if int_contig not in select_tigs:
                select_tigs[int_contig] = entries_dict[int_contig]["len"]
    
    return select_tigs, sum_db_size, x


def plot_ints_per_chrom(Nx_entr_dict, int_dict, x):
    ilog.vlprint("Plotting ints per chromosome", 5)

    def rotation_gage(_int_dict):
        len_ticks = 0
        for i in _int_dict:
            len_ticks += len(i)
        if len_ticks < 60:
            rot = 0
            horiz_aln = "center"
            font_size = 10
        elif len_ticks < 450:
            rot = 45 * math.sin(len_ticks/(90*math.pi))
            horiz_aln = "right"
            if len(_int_dict) < 26:
                font_size = 10
            else:
                font_size = 10*(25/len(_int_dict))
        else:
            rot = 45
            horiz_aln = "right"
            if len(_int_dict) < 26:
                font_size = 10
            else:
                font_size = 10*(25/len(_int_dict))
        return rot, horiz_aln, font_size

    def plot_only_int(_int_dict):
        #plot showing only chromosomes/contigs harboring integrations
        plt.bar(range(len(int_dict)), int_dict.values(), align='center')
        rot_gut = rotation_gage(_int_dict)[0] 
        horiz_aln = rotation_gage(_int_dict)[1]
        fs = rotation_gage(_int_dict)[2]
        plt.xticks(
            range(len(int_dict)),
            int_dict.keys(),
            rotation=rot_gut,
            ha=horiz_aln,
            fontsize=fs
            )
        plt.ylabel(
            "Number of integrations",
            fontdict={"weight": "normal", "size": 14}
            )
        # plt.title("Figure 2: Integrations per contig", loc="left",
        #           fontdict={"weight": "bold", "size": 14},)
        plt.tight_layout()
        plt.savefig("ints_per_chrom.png")
        plt.close()

    def plot_all(Nx_entr_dict, _int_dict):
        #plot showing chromosomes/contigs accounting for >= Nx (98%) of genome
        all_chr = _int_dict
        for cont in Nx_entr_dict:
            if cont in all_chr:
                pass
            else:
                all_chr[cont] = 0
        keys = sorted(Nx_entr_dict.keys(), key=lambda k: Nx_entr_dict[k])
        keys.reverse()
        vals = [all_chr[k] for k in keys]
        plt.bar(range(len(all_chr)), vals, align='center')
        # plt.bar(range(len(vals)), vals, align='center')
        rot_gut = rotation_gage(all_chr)[0]
        horiz_aln = rotation_gage(all_chr)[1]
        fs = rotation_gage(_int_dict)[2]
        plt.xticks(range(
            len(all_chr)), keys, rotation=rot_gut, ha=horiz_aln, fontsize=fs
            )
        plt.ylabel(
            "Number of integrations", fontdict={"weight": "normal", "size": 14}
            )
        # plt.title(f"Figure 3: Integrations per contig (all N{x})", loc="left",
        #           fontdict={"weight": "bold", "size": 14},)
        plt.tight_layout()
        plt.savefig("ints_per_chrom_incl0.png")
        plt.close()

    # plot_all(Nx_entr_dict, int_dict)
    plot_only_int(int_dict)
    plot_all(Nx_entr_dict, int_dict)


def draw_ideogram(
    Nx_entr_dict,
    sum_db_size,
    pot_int_site_objs,
    stroke=20,
    fill_col= "BDBDBD", # "E0E0E0"
    font=20,
    width=840,
    height=600,
    margins=20,
    prefix=""
    ):
    ilog.vlprint("Drawing generic ideogram", 5)
    max_len = width
    graph_contigs = {}
    for contig in Nx_entr_dict:
        if sum_db_size > 0:
            cont_prop = Nx_entr_dict[contig]/sum_db_size
        else:
            cont_prop = 0
        graph_cont_len = max_len * cont_prop
        graph_contigs[contig] = [graph_cont_len]
    # longest_cont = max(graph_contigs.values())
    if len(graph_contigs) > 0:
        adj_stroke = stroke * (stroke/len(graph_contigs))
        adj_font = font * (font/len(graph_contigs))
    else:
        adj_stroke, adj_font = 0, 0
    # if adj_stroke > 100:
    #     step_size = 100
    if adj_font > 42:
        adj_font = 42

    svg = f'{prefix}relative_integration_locations.svg'
    with open(svg, 'w') as pic:
        pic.write('<?xml version = "1.0" ?>\n')
        pic.write(
            f'<svg width="{width}" height="{height}" dpi="72" version="1.1" '
            f'xmlns="http://www.w3.org/2000/svg" '
            f'xmlns:xlink="http://www.w3.org/1999/xlink">\n'
            )
        pic.write(
            f'<rect x="0" y="0" fill="#FFFFFF" stroke="#000000" stroke-width="3"'
            f' width="{width}" height="{height}"/>\n'
            )
        cont_count = 0
        if len(graph_contigs) > 0:
            step_size = (height*0.9/len(graph_contigs)) + adj_stroke/(2*(len(graph_contigs)/1))
        else:
            step_size = 0
        if step_size > 125:
            step_size = 125
        for contig in graph_contigs:
            cont_count += 1
            step_incr = cont_count*step_size
            length = graph_contigs[contig][0]*3.5
            pic.write(
                f'<text x="{length+margins+20}" y="{step_incr + adj_stroke/2}"'
                f' font-size="{adj_font}" >{contig}</text>\n'
                )
            pic.write(
                f'<line stroke="#{fill_col}" stroke-width="{adj_stroke}"'
                f' x1="{margins}" y1="{step_incr}" x2="{length+margins}"'
                f' y2="{step_incr}" stroke-linecap="round"/>\n'
                )
            for site in pot_int_site_objs:
                if contig == site.chrom:
                    # if type(site.approx_loc) == tuple:
                    #     print("site.approx_loc", site.approx_loc)
                    #     ilog.vlprint(
                    #         f"INFO: site.approx_loc of {site.name} is of type tuple", 5
                    #         )
                    #     # site.approx_loc = int(sum(site.approx_loc)/len(site.approx_loc))
                    #     site.approx_loc = site.approx_loc[1]
                    rel_pos = (site.approx_loc/Nx_entr_dict[contig])*length
                    pic.write(
                        f'<line fill="none" stroke="#ff0000" stroke-width="{2}"'
                        f' x1="{rel_pos+margins}" y1="{step_incr-adj_stroke/2}"'
                        f' x2="{rel_pos+margins}" y2="{step_incr+adj_stroke/2}"/>\n'
                        )
        pic.write('</svg>\n')


def draw_generic_ideogram(int_dict, pot_int_objs, args, prefix=""):
    if not any([args.intra, args.polyclonal, args.quick_mm2]):
        try:
            all_tigs = get_genome_info_blast(args.genome)
            rel_tigs = select_contigs(all_tigs, int_dict, x=args.Nx, n=args.ntigs)
        except (KeyError, TypeError) as e:
            ilog.vlprint(e, 10)
            ilog.vlprint(
                "BLAST based retrieval of contig data from database failed. This "
                "is a minor hiccup. Scanning genome reference to retrieve contig "
                "names without BLAST.", 10
                )
            all_tigs = scan_ref_genome(args.genome, int_dict)
            rel_tigs = select_contigs(all_tigs, int_dict, x=args.Nx, n=args.ntigs)
    else:
        all_tigs = scan_ref_genome(args.genome, int_dict)
        rel_tigs = select_contigs(all_tigs, int_dict, x=args.Nx, n=args.ntigs)
    plot_ints_per_chrom(rel_tigs[0], int_dict, rel_tigs[2])
    draw_ideogram(rel_tigs[0], rel_tigs[1], pot_int_objs, prefix=prefix)

