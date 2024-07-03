# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as hi_cl
import random
from il_logging import Ilogger

ilog = Ilogger()
ilog.module = __name__


def blast_reads_all_vs_all(args):
    """All-vs-all alignment of candidate reads using BLAST"""
    ilog.vlprint("Performing all-vs-all alignment for read clustering", 3)
    cr_wo_ints = "cand_reads_without_integration.fasta"
    bdb_fext = ['.ndb', '.nal', '.nsq', '.00.nsq']
    if any([os.path.exists(cr_wo_ints + ext) for ext in bdb_fext]):
        ilog.vlprint('Read-DB exists. Skipping generation of read-db.', 2)
    elif os.path.getsize(cr_wo_ints)==0:
        ilog.vlprint('No candidate reads. Skipping generation of read-db.', 2)
    else:
        ilog.vlprint('Generating BLAST-DB from reads file', 1)
        subprocess.run(
                    [
                    args.makeblastdb,
                    "-in", cr_wo_ints,
                    "-dbtype", "nucl",
                    "-parse_seqids",
                ],
                check=True, capture_output=True
            )
    if os.path.getsize(cr_wo_ints)>0:  
        subprocess.run(
                [
                args.blastn,
                "-db", cr_wo_ints,
                "-query", cr_wo_ints,
                "-out", "ava_res.csv",
                "-outfmt", "10 qaccver qlen saccver slen length qstart qend sstart \
                send pident mismatch gapopen evalue bitscore qcovhsp qcovus",
                "-num_threads", "2"
                ],
                check=True, capture_output=True
                )
    else:
        ilog.vlprint('No candidate reads. Skipping ava-blast.', 2)
        open("ava_res.csv", "w").close()
                     

def gen_ava_al_dict():
    """Generate a dictionary with query_reads as primary keys and dicts as 
       values with target sequences as keys and the query_to_target alignment-
       results as value for clustering"""
    ilog.vlprint("Generating merged_als dict for clustering", 3)
    als = []
    with open("ava_res.csv", "r") as ava:
        for line in ava:
            sl = line.split(",")
            str_ids = [0, 2]
            for i in range(len(sl)):
                if i not in str_ids:
                    sl[i] = float(sl[i])
            als.append(sl)

    # merged_als is a dict, with query_reads as primary keys and dicts as vals,
    # which in turn have subj/target_sequences as keys and the respective 
    # query_to_target alignment-results as value
    merged_als = {}
    for line in als:
        if line[0] not in merged_als:
            merged_als[line[0]] = {line[2]: line[3:]}
            for i in [1, 2, 3, 4, 5, 6, 9]:
                merged_als[line[0]][line[2]][i] = list([line[3:][i]])
        else:
            # for al in merged_als:
            if line[2] not in merged_als[line[0]]:
                merged_als[line[0]][line[2]] = line[3:]
                for i in [1, 2, 3, 4, 5, 6, 9]:
                    merged_als[line[0]][line[2]][i] = list([line[3:][i]])
            else:
                for i in [1, 2, 3, 4, 5, 6, 9]:
                    merged_als[line[0]][line[2]][i].append(line[3:][i])
                merged_als[line[0]][line[2]][7] += line[3:][7]
                merged_als[line[0]][line[2]][8] += line[3:][8]
                merged_als[line[0]][line[2]][10] += line[3:][10]
                merged_als[line[0]][line[2]][11] += line[3:][11]

    return merged_als


def filter_merged_for_clustering(_merged_als, thr_bitscore=200):
    """Filter merged alignments based on cumulative bitscore"""
    ilog.vlprint("Filtering merged_als dict", 3)
    filt_merged = {}
    for key in _merged_als:
        filt_merged[key] = {}
        for sub_key in _merged_als[key]:
            if _merged_als[key][sub_key][-3] > thr_bitscore:
                filt_merged[key][sub_key] = _merged_als[key][sub_key]

    return filt_merged


class AvaTarget:
    def __init__(self, query, target, tlen, al_len, qstart, qend, tstart, tend,
                 pident, mismatch, gapopen, evalue, bitscore, qcovhsp, qcovus):
        self.query = query
        self.target = target
        self.tlen = tlen
        self.al_len = al_len
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend
        self.pident = pident
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.evalue = evalue
        self.bitscore = bitscore
        self.qcovhsp = qcovhsp
        self.qcovus = qcovus


class AvaQuery:
    def __init__(self, query, targets=[]):
        self.query = query
        self.targets = targets
        self.bitscore_best = None
        self.read_names = set()
        self.read_names_plus = None
        self.qlen = None

    def __str__(self):
        return f"AvA-Query_{self.query}"

    def __repr__(self):
        return f"AvA-Query_{self.query}"

    def set_qlen(self):
        for target in self.targets:
            if target.target == self.query:
                self.qlen = float(target.tlen)

    def best_match(self):
        """Find best match within target sequences based on bitscore"""
        bitscore_matches = []
        for target in self.targets:
            match_tup = (target.target, target.bitscore)
            bitscore_matches.append(match_tup)
        if len(bitscore_matches) > 1:
            nonself_best = sorted(bitscore_matches, key=lambda x: x[1])[-2]
        else:
            nonself_best = False
        self.bitscore_best = nonself_best

    def calc_cum_bitscore(self, selected_targets=[]):
        """Calculate the cumulative bitscore for the intersection of target-seqs
         between self and other query"""
        cum_bitscore = 0
        for target in self.targets:
            if target.target in selected_targets:
                cum_bitscore += target.bitscore

        return cum_bitscore

    def populate_read_names(self):
        self.read_names.add(self.query)
        for target in self.targets:
            self.read_names.add(target.target)
        self.read_names_plus = self.read_names
        self.read_names_plus.add(self.query)


def gen_target_and_query_instances(_filt_merged):
    """Generate target and query object instances"""
    ilog.vlprint("Generating target and query instances", 3)
    nested_al_obj = []
    for qname in _filt_merged:
        target_lst = []
        for tname in _filt_merged[qname]:
            tlist = _filt_merged[qname][tname]
            tobj = AvaTarget(
                qname,
                tname,
                tlist[0],
                tlist[1],
                tlist[2],
                tlist[3],
                tlist[4],
                tlist[5],
                tlist[6],
                tlist[7],
                tlist[8],
                tlist[9],
                tlist[10],
                tlist[11],
                tlist[12]
                )
            target_lst.append(tobj)
        aqobj = AvaQuery(qname, target_lst)
        nested_al_obj.append(aqobj)
    for qname in nested_al_obj:
        qname.best_match()
        qname.populate_read_names()
        qname.set_qlen()

    return nested_al_obj

class PotAvaCluster:
    def __init__(self, name, ava_queries=[]):
        self.name = name
        self.ava_queries = ava_queries
        self.read_names = set()
        self.cum_bitscores = None
        self.overlaps = {}
        self.parents = False
        self.fuse_to = False
        self.fusable = True
        self.block_fusion = []
        self.multiway_merger = False
        self.rm_merged = False

    def merge(self, obj2):
            self.ava_queries.extend(obj2.ava_queries)

            try:
                self.cum_bitscores = dict(
                    **self.cum_bitscores, **obj2.cum_bitscores
                    )
                try:
                    self.overlaps = dict(**self.overlaps, **obj2.overlaps)
                except TypeError:
                    for group in obj2.overlaps:
                        # may be removed providing additional testing; 
                        for query_read in obj2.overlaps[group]:
                            if group in self.overlaps:
                                self.overlaps[group].update(obj2.overlaps[group])
                            else:
                                self.overlaps[group] = obj2.overlaps[group]
                self.read_names = self.read_names.union(obj2.read_names)
                self.parents = [self.name, obj2.name]
                self.name = self.name + "_" + obj2.name
                self.fuse_to = False
                obj2.fuse_to = False
                obj2.fusable = False
                obj2.rm_merged = True
            except TypeError as e:
                print("Attempted merging of clusters failed")
                print(e)

    def populate_read_names(self):
        for query in self.ava_queries:
            self.read_names.add(query.query)

    def get_ovl_other_groups(self):
        for query_name in self.cum_bitscores:
            for tup in self.cum_bitscores[query_name]:
                if not tup[1]:
                    if tup[2] > 0:
                        if tup[0] not in self.overlaps:
                            self.overlaps[tup[0]] = {
                                query_name: dict(cum_bs=tup[2], av_bs=tup[3])
                                }
                        else:
                            self.overlaps[tup[0]][query_name] = dict(
                                cum_bs=tup[2], av_bs=tup[3]
                                )

    def sum_cbs(self, other_group):
        """Calculate sum of the cumulative bitscores (cumcum-bs) of all queries 
        overlapping with other group"""
        all_cbs = []
        for read in self.overlaps[other_group]:
            all_cbs.append(self.overlaps[other_group][read]["cum_bs"])
        return sum(all_cbs)

    def sum_avbs(self, other_group):
        """Calculate sum of the average bitscores (cumav-bs) of all queries 
        overlapping with other group"""
        all_avbs = []
        for read in self.overlaps[other_group]:
            all_avbs.append(self.overlaps[other_group][read]["av_bs"])
        return sum(all_avbs)

    def cull_ovl(self, bs_ovl_thr=500):
        # remove other groups that have less than bs_ovl_thr total cumulative 
        # bitscore overlap with
        self.overlaps = {
            k: v for (k, v) in self.overlaps.items() if self.sum_cbs(k) > bs_ovl_thr
            }
        # set fuse_to flag to group with highest cumulative and average bitscore
        #  values if overlap exists
        if len(self.overlaps) > 0:
            # for this first create "inverted" dict with cumulated cum_cbs 
            # values as keys and group_names as values
            inv_cumbs_groups = {
                self.sum_cbs(k): k for (k, v) in self.overlaps.items()
                }
            # then set fuse_to flag to the value of the biggest key in the inv-dict
            self.fuse_to = inv_cumbs_groups[max(inv_cumbs_groups)]
            # may become optional: if the summed average bit_score values are 
            # below thr.(500 for now) cancel fuse-flag
            if self.sum_avbs(self.fuse_to) < bs_ovl_thr:
                self.fuse_to = False
            # create a list of prohibited fusion partner groups based on low 
            # average cum_bitscore
            self.block_fusion = [
                k for k in self.overlaps if self.sum_avbs(k) < (self.sum_cbs(k)/20)
                ]
            if len(self.overlaps) > 1:
                # use set comprehension to include all group names from overlaps
                self.multiway_merger = {s for s in self.overlaps}
                # add self.name to be able to compare identity between sets with
                #  other groups
                self.multiway_merger.add(self.name)
    
    def report_self(self):
        print(self.name)
        print("Number of clustered reads: ", len(self.read_names))


def gen_read_cluster(ava_query_objs):
    """Generate read clusters from all-vs-all query objects"""
    ilog.vlprint("Generating read cluster", 3)
    # remove queries that have no matches (other than self), as those can't be 
    # clustered
    ava_query_objs = [query for query in ava_query_objs if len(query.targets) > 1]
    # make_adjustable: remove all queries with length < 2500 bp
    ava_query_objs = [query for query in ava_query_objs if query.qlen > 2500]
    # experimental: add a random component to check how initial nucleation of 
    # groups affects later cluster growth consider using a network/ interconnect-edness
    #  parameter, which indicates with how many other reads a subject read is 
    # related; such a parameter could prove a good guidance for selecting the 
    # optimal read for initiation of cluster growth (at least better then 
    # merely length) and may make random component obsolete
    random.shuffle(ava_query_objs)

    def bbb_group_queries(ava_query_objs, qbins):
        """Add queries to qbins; create new bin if best_match not in qbins, else
         add to bin with best match"""
        while len(ava_query_objs) > 0:
            best_match = ava_query_objs[0].bitscore_best[0]
            in_qbins = False
            for qbin in qbins:
                for query in qbin:
                    if query.query == best_match:
                        qbin.append(ava_query_objs[0])
                        ava_query_objs.pop(0)
                        in_qbins = True
            if not in_qbins:
                qbin = list()
                qbin.append(ava_query_objs[0])
                qbins.append(qbin)
                ava_query_objs.pop(0)
        return qbins

    def extract_unique(qbins):
        """Move single bins back to reservoire for next iteration of 
        bbb_group_queries"""
        single_qobjs = [qbin[0] for qbin in qbins if len(qbin) == 1]
        qbins = [qbin for qbin in qbins if len(qbin) > 1]
        return qbins, single_qobjs

    def iterative_grouping(ava_query_objs):
        """run iterations of bbb_group_queries until all queries are assigned to
         a cluster"""
        qbins = []
        ic = 0
        while len(ava_query_objs) > 0 and ic < 3:
            ic += 1
            qbins = bbb_group_queries(ava_query_objs, qbins)
            qbins, ava_query_objs = extract_unique(qbins)
        ilog.vlprint(
            f"Ran {ic} iterations of read grouping to form potential clusters", 1
            )
        ilog.vlprint(
            f"Number of potential Clusters: {len(qbins)}", 6, logging=False
            )
        return qbins

    pclust = iterative_grouping(ava_query_objs)
    return pclust


def gen_potclust_instances(pclust):
    """Generate instances of potential read clusters resulting from all-vs-all 
    alignment"""
    ilog.vlprint("Generating potential cluster instances", 3)
    pcc = 0
    pci_lst = []
    for pc in pclust:
        pcc += 1
        pc_name = "pC" + str(pcc)
        pci = PotAvaCluster(pc_name, pc)
        pci.populate_read_names()
        pci_lst.append(pci)
    return pci_lst


def chk_cum_bitscore_ovl(pci_lst, devdeb=False):
    """Evaluate if the query_read is grouped together with reads giving maximum 
    cumulative bitscore and if significant similarity exists to other read group
     indicating groups that have to be merged"""
    ilog.vlprint("Checking cumulative bitscores", 3)
    occurs = {}

    # for each potentialCLuster(read-group) instance
    for pci in pci_lst:
        cum_bitscores = {}
        # for every query in this instance
        for query in pci.ava_queries:
            # count how often this query occurs in all groups
            if query.query not in occurs:
                occurs[query.query] = 1
            else:
                occurs[query.query] += 1
            # check in pci-list fin every potentialCLuster(read-group) instance
            for group in pci_lst:
                # if it is the cluster/group the query is contained in mark as 
                # True
                if query.query in group.read_names:
                    isin_group = True
                # if it is any other cluster/group than the query is assigned to
                # mark as False
                else:
                    isin_group = False
                # collect names of all reads contained in the other pC/group
                reads_in_group = list(group.read_names)
                # calculate the bit-score overlap of the query with all those 
                # reads (based on the targets in query)
                cum_score = query.calc_cum_bitscore(reads_in_group)
                av_cum_score = cum_score/len(reads_in_group)
                if query.query not in cum_bitscores:
                    cum_bitscores[query.query] = [
                        (group.name, isin_group, cum_score, av_cum_score)
                        ]
                else:
                    cum_bitscores[query.query].append(
                        (group.name, isin_group, cum_score, av_cum_score)
                        )
        pci.cum_bitscores = cum_bitscores
        pci.get_ovl_other_groups()
        pci.cull_ovl()

    if devdeb:
        for pci in pci_lst:
            if len(pci.overlaps) > 0:
                ilog.vlprint(f"Culled cbs {pci.name}:", 2, logging=False)
                ilog.vlprint(pci.overlaps, 2, logging=False)

        for pci in pci_lst:
            ilog.vlprint(
                f"Cumulative bitscore infos for group {pci.name}:",
                2, logging=False
                )
            for query_name in pci.cum_bitscores:
                assigned_cluster = ""
                cbs_own_cl = int()
                for tup in pci.cum_bitscores[query_name]:
                    if not tup[1]:
                        if tup[2] > 0:
                            ilog.vlprint(
                                f"{query_name} has overlap with other cluster "
                                f"{tup[0]} ;cum_ovl = {tup[2]}, av_ovl = {tup[3]}",
                                2, logging=False)
                    else:
                        assigned_cluster = tup[0]
                        cbs_own_cl = tup[2]
                        av = tup[3]
                print(
                    f"{query} is assigned to cluster {assigned_cluster} with a "
                    f"cumulative bit_score of {cbs_own_cl} and an average "
                    f"bit_score of {av}"
                    )

    return pci_lst


def merge_ovl_groups(obj_lst, min_read_supp):
    """merge groups with the highest scoring overlapping group"""
    ilog.vlprint("Merging overlapping groups", 3)
    # perform multiway-merger between more than two groups
    for obj in obj_lst:
        mwm_temp = []
        if obj.multiway_merger:
            for any_obj in obj_lst:
                if any_obj.multiway_merger == obj.multiway_merger:
                    mwm_temp.append(any_obj)
            if len(mwm_temp) == len(obj.multiway_merger):
                for index in range(len(mwm_temp))[1:]:
                    ilog.vlprint(
                        f"Merge {mwm_temp[0].name} to {mwm_temp[index].name}",
                        6, logging=False
                        )
                    mwm_temp[0].merge(mwm_temp[index])
                    ilog.vlprint(
                        f"Merged Result = {mwm_temp[0].name} "
                        f"Parents: {mwm_temp[0].parents}", 6, logging=False
                        )
                    mwm_temp[index].multiway_merger = False
                    mwm_temp[index].merge_to = False
                    mwm_temp[index].fusable = False
                mwm_temp[0].multiway_merger = False
                mwm_temp[0].merge_to = False
                mwm_temp[0].fusable = False

    # perform simple merger between two-groups based on the highest scoring 
    # overlapping group
    for obj in obj_lst:
        if obj.fusable and obj.fuse_to:
            for other in obj_lst:
                if other.name == obj.fuse_to and other.fusable:
                    if obj.name not in other.block_fusion:
                        obj.merge(other)

    obj_lst = [obj for obj in obj_lst if not obj.rm_merged]

    # rename merged groups with nomenclature
    cc = 0
    for cluster in obj_lst:
        cc += 1
        cluster.name = "Cl" + str(cc)
        ilog.vlprint(
            f"{cluster.name} parents: {cluster.parents}", 4, logging=False
            )

    # filter for min_read_support
    obj_lst = [
        cluster for cluster in obj_lst if len(cluster.read_names) >= min_read_supp
        ]

    return obj_lst


def write_ava_cluster_csv(obj_lst):
    ilog.vlprint("Writing ava_read_cluster.csv", 3)
    with open("ava_read_cluster.csv", "w") as out_csv:
        out_csv.write("cluster, #reads, read names\n")
        for cluster in obj_lst:
            out_csv.write(
                ",".join([cluster.name, str(len(cluster.read_names)),
                ";".join(cluster.read_names)]) + "\n"
                )

def dump_reads_by_cluster(obj_lst):
    ilog.vlprint(f"A total number of {len(obj_lst)} read_cluster were formed.", 1)
    clust_names = [obj.name for obj in obj_lst]
    clust_read_dict = {
        ">" + readname + "\n": obj.name for obj in obj_lst for readname in obj.read_names
        }
    outfiles = dict()
    for cname in clust_names:
        outfile = open(cname + "_reads.fasta", "w")
        outfiles[cname] = outfile
    outfiles["not_assignable"] =  open("Cluster_not_assignable_reads.fasta", "w")
    
    
    with open("Integration_candidate_reads.fasta", "r") as int_reads:
        assigned = 0
        unassigned = 0
        total_reads = 0
        for line in int_reads:
            if line.startswith(">"):
                total_reads += 1
                try:
                    cluster = clust_read_dict[line]
                    assigned += 1
                except:
                    cluster = "not_assignable"
                    unassigned += 1
                seq = next(int_reads)
                out = outfiles[cluster]
                out.write(line)
                out.write(seq)

    for outfile in outfiles:
        outfiles[outfile].close() 


def associate_cluster_with_potintsite(clust_objs, int_sites):
    for clust in clust_objs:
        clust.intsite_links = {site.name: 0 for site in int_sites}
        for read in clust.read_names:
            for site in int_sites:
                if read in site.read_names:
                    clust.intsite_links[site.name] += 1
        clust.intsite_links = {
            k: v for k,v in clust.intsite_links.items() if v > 0
            }
        ilog.vlprint(f"Linkages cluster {clust.name} {clust.intsite_links}", 2)
    return clust_objs


def heatmap_ava_clust(obj_lst, _merged_als):
    ilog.vlprint("Generating custom heatmap", 3)
    x_ks = []
    x_ls = []
    for cl in obj_lst:
        for read in cl.read_names:
            x_ks.append(read)
            x_ls.append(cl.name + "_" + read[:8])
    x_keys, y_keys = tuple(x_ks), tuple(x_ks)
    x_labels, y_labels = tuple(x_ls), tuple(x_ls)

    # get query sequence lengths
    qlens = {}
    for query in _merged_als:
        for subject in _merged_als[query]:
            if query == subject:
                qlens[query] = _merged_als[query][subject][0]

    matrix = np.zeros((len(x_keys), len(y_keys)))
    x = -1
    for xkey in x_keys:
        x += 1
        y = -1
        for ykey in y_keys:
            y += 1
            try:
                # Use for plotting of qcovus
                matrix[x][y] = _merged_als[xkey][ykey][-1]
                ################################################################
                ##Alternative plotting options##
                # Use for plotting of qcovhsp
                # matrix[x][y] = _merged_als[xkey][ykey][-2]
                # Use for plotting of bitscore
                # matrix[x][y] = _merged_als[xkey][ykey][-3]
                # Use for plotting of bitscore/qlen ratio
                # matrix[x][y] = _merged_als[xkey][ykey][-3]/(qlens[xkey]*2)
                # Use for plotting of qcovus-bitscore compound value
                # matrix[x][y] = (_merged_als[xkey][ykey][-1] + (
                # _merged_als[xkey][ykey][-3]/(qlens[xkey]*2))*100
                # )/2
                ################################################################

            except KeyError:
                matrix[x][y] = 0
                continue

    fig, ax = plt.subplots()
    im = ax.imshow(matrix, cmap="bwr")

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    # for plotting of qcovus
    cbar.ax.set_ylabel("qcovus [%]", rotation=-90, va="bottom")
    ##########################################################################
    ##Alternative plotting options##
    # Use for plotting of qcovhsp
    # cbar.ax.set_ylabel("qcovhsp [%]", rotation=-90, va="bottom")
    # #for plotting of bitscore
    # cbar.ax.set_ylabel("bitscore", rotation=-90, va="bottom")
    # #for plotting of bitscore/qlen ratio
    # cbar.ax.set_ylabel("bitscore/qlen*2", rotation=-90, va="bottom")
    # Use for plotting of qcovus-bitscore compound value
    # cbar.ax.set_ylabel("qcovus_normbitscore av.", rotation=-90, va="bottom")
    ##########################################################################

    # show all ticks
    ax.set_xticks(np.arange(len(x_keys)))
    ax.set_yticks(np.arange(len(y_keys)))
    # label ticks with respective names
    ax.set_xticklabels(x_labels, fontdict={'fontsize': 2})
    ax.set_yticklabels(y_labels, fontdict={'fontsize': 2})
    # Rotate tick labels and set alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    ax.set_title(
        "Query Coverage Per Unique Subject between integration candidate reads"
        )

    fig.tight_layout()
    plt.savefig("read_cluster_il-ava-clust.svg")
    plt.savefig("read_cluster_il-ava-clust.png")
    plt.close()


def heatmap_scipy_clust(obj_lst, _merged_als):
    ilog.vlprint("Generating scipy heatmap", 3)
    x_ks = []
    for cl in obj_lst:
        for read in cl.read_names:
            x_ks.append(read)
    x_keys, y_keys = tuple(x_ks), tuple(x_ks)

    matrix = np.zeros((len(x_keys), len(y_keys)))
    x = -1
    for xkey in x_keys:
        x += 1
        y = -1
        for ykey in y_keys:
            y += 1
            try:
                # Use for plotting of qcovus
                matrix[x][y] = _merged_als[xkey][ykey][-1]
            except KeyError:
                matrix[x][y] = 0
                continue

    # Dendrogram that comes to the left
    fig = plt.figure(figsize=(15, 15))
    ax1 = fig.add_axes([0.1, 0.1, 0.2, 0.6])
    link = hi_cl.linkage(matrix, method='single') # centroid
    tree = hi_cl.dendrogram(link, orientation='left')
    for k in ax1.spines:
        ax1.spines[k].set_visible(False)
    ax1.set_xticks([])
    ax1.set_yticks([])

    # top side dendrogram
    ax2 = fig.add_axes([0.3, 0.70, 0.6, 0.2])
    link = hi_cl.linkage(matrix, method='single')
    tree = hi_cl.dendrogram(link)
    for k in ax2.spines:
        ax2.spines[k].set_visible(False)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # main heat-map
    axmatrix = fig.add_axes(
        [0.3, 0.1, 0.6, 0.6],
        yticklabels=[],
        yticks=[],
        xticklabels=[],
        xticks=[]
        )
    idx = tree['leaves']
    # use numpy indexing/slicing to first sort rows according to the order of 
    # indices in idx1
    matrix = matrix[idx, :]
    # use numpy indexing/slicing to then sort columns according to the order of 
    # indices in idx2
    matrix = matrix[:, idx]
    # show heat-map image (use imshow; matshow will re-introduce xticks)
    im = axmatrix.imshow(matrix, aspect='auto', origin='lower', cmap="bwr")

    plt.savefig("read_cluster_ava_scipy.svg")
    plt.savefig("read_cluster_ava_scipy.png")
    plt.close()


def ava_cluster(min_read_supp, args):
    """Run all-vs-all alignment and clustering workflow"""
    ilog.vlprint("Running all-vs-all alignment and clustering workflow", 3)
    if min_read_supp < 3:
        ilog.vlprint(
            "Read clustering cannot be performed with minimum read support "
            "lower than 3. Skipping clustering.", 0 
            )
    else:
        blast_reads_all_vs_all(args)
        merged_als = gen_ava_al_dict()
        filt_merged = filter_merged_for_clustering(merged_als)
        ava_al_objs = gen_target_and_query_instances(filt_merged)
        pclust = gen_read_cluster(ava_al_objs)
        pci_lst = gen_potclust_instances(pclust)
        ovl_groups = chk_cum_bitscore_ovl(pci_lst)
        merge_groups = merge_ovl_groups(ovl_groups, min_read_supp)
        ovl_groups = chk_cum_bitscore_ovl(merge_groups)
        merge_groups = merge_ovl_groups(ovl_groups, min_read_supp)
        ovl_groups = chk_cum_bitscore_ovl(merge_groups)
        merge_groups = merge_ovl_groups(ovl_groups, min_read_supp)
        ilog.vlprint(
            f"Integration candidate reads clustered into {len(merge_groups)}"
            f" groups", 0
            )
        write_ava_cluster_csv(merge_groups)
        dump_reads_by_cluster(merge_groups)
        heatmap_ava_clust(merge_groups, merged_als)
        if len(merge_groups)>0:
            heatmap_scipy_clust(merge_groups, merged_als)
        else:
            open("read_cluster_ava_scipy.svg", "w").close()
            open("read_cluster_ava_scipy.png", "w").close()
        return merge_groups


if __name__ == "__main__":
    # The following two function calls are only required if ran independently of 
    # intloc main script
    # cand_reads = gen_cand_read_instances("candidate_reads.fbr")
    # write_reads_wo_integration(cand_reads)
    from integration_locator import BlastCandidateRead
    from integration_locator import gen_cand_read_instances, report_time
    report_time("All-vs-all read clustering")
    ava_cluster(1)
    report_time("All-vs-all read clustering", end=True)
