# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.


import os
import sys
from time import strftime, localtime, perf_counter
from collections import defaultdict
import random
import argparse
import statistics as stat


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Integrator will incorporate a supplied DNA sequence (integ"
                    "ration) a specified number of times into an acceptor seque"
                    "nce (or collection of sequences, typically a reference gen"
                    "ome) at either randomly chosen locations or in a targeted "
                    "manner thereby generating a copy of the acceptor sequence"
                    " carrying these integrations. This can be used to simulate"
                    " reads e.g. for testing of intloc.", 
                    add_help=True
                    )

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument(
        "-i", "--integration", required=True, 
        help="single-line fasta file with integrating sequence"
                    )
    required_args.add_argument(
        "-a", "--acceptor", required=True,
        help="fasta file of acceptor sequence (e.g. ref.genome)"
                    )
    parser.add_argument(
        "-n", "--number", default=1, help="number of times the integration will"
        " be pasted into the acceptor sequence"
                        )
    parser.add_argument(
        "--targeted", default=False, help="provide path to a semicolon"
        " separated file containing one location site per line described by:"
        " chromosome_header;base_position_after_which_to_insert"
                        )
    args = parser.parse_args()
    return args


def parse_acceptor_seq(acceptor):
    """Scan assembly into which the integrations will be inserted"""
    contigs = dict()
    curr_cont = ""
    n_stretch = False
    curr_ns_start = None
    n_regions = defaultdict(list)
    with open(acceptor, "r") as acc:
        for line in acc:
            if line.startswith(">"):
                if n_stretch:
                    n_stretch = False
                    cns_end = contigs[curr_cont]
                    n_regions[curr_cont].append((curr_ns_start, cns_end))
                curr_cont = line.strip()
                contigs[curr_cont] = 0
            else:
                if n_stretch:
                    bases = ["A","T","G","C","a","t","g","c"]
                    chk_end = [base for base in line if base in bases]
                    if any(chk_end):
                        n_stretch = False
                        cns_end = min([line.index(base) for base in chk_end])
                        cns_end = cns_end + contigs[curr_cont]
                        n_regions[curr_cont].append((curr_ns_start, cns_end))
                else:
                    if "N" in line:
                        n_stretch = True
                        curr_ns_start = contigs[curr_cont] + line.index("N")
                contigs[curr_cont] += len(line.strip())
        if n_stretch:
                    n_stretch = False
                    cns_end = contigs[curr_cont]
                    n_regions[curr_cont].append((curr_ns_start, cns_end))

    stats = [contigs[tig] for tig in contigs]
    total_bp = sum(stats)
    gaps = [reg for reg_lst in n_regions.values() for reg in reg_lst]
    gap_lens = [reg[1]-reg[0] for reg in gaps]
    total_ns = sum([reg[1]-reg[0] for reg in gaps])

    print("Acceptor Seq. statistics:")
    print("No. of contigs", len(stats))
    print("Average contig length:", stat.mean(stats))
    print("longest contig", max(stats))
    print("Total base count", total_bp)
    print("Number of gaps (N-stretches) in genome", len(gaps))
    print("Total N count in genome", total_ns)
    if len(gap_lens) > 0:
        print("Average gap length", sum(gap_lens)/len(gap_lens))

    return contigs, total_bp, n_regions


def pick_int_sites(num, contigs, total_bp, out_name, n_regions):
    """Randomly pick genomic sites where integrations will later be inserted"""
    
    def random_pick_sites(num, contigs, total_bp):
        "randomly select genomic sites"
        int_sites = []
        for n in range(num):
            conts = list(contigs.keys())
            weights = [contigs[tig]/total_bp for tig in conts]
            target_tig = str(*random.choices(conts, weights=tuple(weights)))
            target_pos = random.choice(range(contigs[target_tig]))
            int_sites.append((n+1, target_tig, target_pos))
        return int_sites

    def chk_for_gap_sites(int_sites, n_regions):
        """Check if any of the picked sites are lying within an assembly gap"""
        print("Checking for ints placed in gaps in assembly")
        in_gaps = list()
        for ins in int_sites:
            id, chrom, loc = ins[0], ins[1], ins[2]
            for gap in n_regions[chrom]:
                if gap[0]-2500 <= loc <= gap[1]+2500:
                    in_gaps.append(id)
        do_repick = len(in_gaps)
        for site in in_gaps:
            print(
                f"Pick {site} will be excluded, because it is in an assembly-gap."
                )
        print(f"For {do_repick} ints a new random pick is required")
        int_sites = [ins for ins in int_sites if ins[0] not in in_gaps]
        return do_repick, int_sites, in_gaps

    int_sites = random_pick_sites(num, contigs, total_bp)
    repick, int_sites, ids = chk_for_gap_sites(int_sites, n_regions)
    new_picks = []
    while repick > 0:
        new = random_pick_sites(repick, contigs, total_bp)
        new = [(ids[new.index(pick)],pick[1],pick[2]) for pick in new]
        repick, new, ids = chk_for_gap_sites(new, n_regions)
        new_picks.extend(new)
    
    for site in new_picks:
        print(f"{site} will be included, to replace rejected gap sites.")
    
    int_sites.extend(new_picks)
    int_sites = sorted(int_sites, key=lambda x: (x[1], x[2]))

    int_dict = defaultdict(list)
    for ins in int_sites:
        if ins[1] in int_dict:
            int_dict[ins[1]].append(ins[2])
        else:
            int_dict[ins[1]] = [ins[2]]
    report = "integrator_report_" + out_name + ".csv"
    with open(report, "w") as rep:
        for ins in int_sites:
            contig = ins[1].replace(",", ";")
            rep.write(f"{ins[0]}, {contig}, {ins[2]}\n")

    return int_dict


def load_target_sites(targets_tsv, contigs, out_name):
    int_dict = defaultdict(list)
    with open(targets_tsv, "r") as targets:
        for line in targets:
                ls = line.strip()
                lss = ls.split(";")
                if len(ls) == 0:
                    pass
                else:
                    if len(lss) == 2:
                        int_dict[lss[0].strip()].append(int(lss[1].strip()))
                    else:
                        print(
                        "Target site file format incorrect. Use semicolon as"
                        " separator. Don't use headers or comments."
                        " Only use two values per line: 1.) The full"
                        " chromosome name, corresponding to the fasta header in"
                        " the genome assembly. 2.) A number indicating the base"
                        " position after which the integration is supposed to" 
                        " be inserted. Separated from each other by a semicolon"
                        ". Example: "
                        " >NC_000021.9 Homo sapiens chromosome 21, GRCh38.p12"
                        " Primary Assembly;1678325"
                        )
        
    # check compatibility of loaded chroms and positions with tigs in acceptor
    for chrom in int_dict:
        if chrom not in contigs:
            print(f"Contig {chrom} is not present in acceptor assembly.")
            print("Available contigs/chromosomes are:")
            for tig in contigs:
                print(tig)
            print("Please edit target sites csv file accordingly.")
            print("Aborting integration_locator.")
            sys.exit()
        for pos in int_dict[chrom]:
            if pos > contigs[chrom]:
                print(f"Position {pos} exceeds length of contig {chrom} with"
                       " length {contigs[chrom]}")
                print("Please edit target sites tsv file accordingly.")
                print("Aborting integration_locator.")
                sys.exit()

    report = "integrator_report_" + out_name + ".csv"
    with open(report, "w") as rep:
        int_ct = 0
        for chrom in int_dict:
            chr_name = chrom.replace(",", ";")
            for pos in int_dict[chrom]:
                int_ct += 1
                rep.write(f"{int_ct}, {chr_name}, {pos}\n")
        
    return int_dict


def insert_ints(int_dict, acceptor, integration, out_name):
    """Insert integrating sequence at previously picked genomic locations"""
    with open(integration) as int_file:
        header_ct = 0
        seq_name = ""
        int_seq = ""
        l_ct = 0
        bases = ["a", "t", "c", "g", "A", "T", "C", "G", "N", "n"]
        for line in int_file:
            l_ct += 1
            if header_ct > 1:
                print(
                    "Fasta file with integrating sequence contains more than" 
                    " one sequence item. Please provide a fasta-file with only"
                    " a single sequence to integrate."
                    )
                exit()
            elif line.startswith(">"):
                header_ct += 1
                seq_name = line.strip()
            elif line[0] in bases:
                int_seq = int_seq + line.strip()
            elif line[0] not in bases:
                if len(line.strip()) > 0:
                    print(
                        f"Unexpected base/character in line {l_ct} in integrati"
                        f"on sequence file. Aborting."
                        )
                    exit()
            else:
                pass

    report = "integrator_report_" + out_name + ".csv"
    with open(report, "a") as rep:
        rep.write("\n")
        if len(seq_name) > 0:
            rep.write(f"Integrating sequence name: {seq_name}\n")
            rep.write(f"Integrating sequence length: {len(int_seq)}")
        else: 
            rep.write(
                "Integrating sequence name: no integration name was provided\n"
                )
            rep.write(f"Integrating sequence length: {len(int_seq)}")

    outfile = open(out_name + ".ints", "w")
    with open(acceptor, "r") as acc:
        lc = 0
        curr_lines = list()
        contigs = dict()
        inserts = 0
        for line in acc:
            lc += 1
            if line.startswith(">"):
                    if lc > 1:
                        contigs[header].seq = curr_lines
                        curr_lines = list() 
                        contigs[header].join_seq_lines()
                        inserts += contigs[header].insert_ints(int_seq)
                        contigs[header].write_out(outfile)
                    header = line.strip()
                    ints = int_dict[header]
                    tig = Contig(header, ints)
                    contigs[header] = tig
            else:
                curr_lines.append(line.strip())
        contigs[header].seq = curr_lines
        contigs[header].join_seq_lines()       
        inserts += contigs[header].insert_ints(int_seq)
        contigs[header].write_out(outfile)
    print(f"A total of {inserts} integrations were inserted into the genome.")
    outfile.close()


class Contig():
    """Class for loading, modifying and writing of contigs in assembly"""
    def __init__(self, header, ints):
        self.header = header
        self.ins = sorted(ints)
        self.mask = self.ins
        self.seq = list()
        self.inserted = ""
        # print(f"Created contig {self.header}")
    
    def join_seq_lines(self):
        self.seq = "".join(self.seq)

    def insert_ints(self, intseq):
        """Insert integrating sequence into locations specified by self.ins"""
        if len(self.ins) > 0:
            prev_int = int()
            for ins in self.ins:
                # in case of the extremly unlikely event that the same location
                # on the same chromosome is chosen multiple times, this will
                # result in direct head-to-tail insertions of the integrating
                # sequence without intervening bases from the acceptor sequence,
                # because [:ins - prev_int] results in [:0] = "".
                self.inserted = self.inserted + self.seq[:ins - prev_int] + intseq
                self.seq = self.seq[ins-prev_int:]
                prev_int = ins
            self.inserted = self.inserted + self.seq
            self.seq = ""
        else:
            self.inserted = self.seq
            self.seq = ""
        print(
            f"succesfully inserted {len(self.ins)} integrations into {self.header}"
            )
        return len(self.ins)
    
    def mask_seq(self, delete=False, logger=False):
        """Mask the parts of the sequence specified by self.mask"""
        # avoid delete option, since it changes length of genome and therefore
        # will confuse all subsequent calculations based on reference position
        bases_masked = 0
        self.inserted = self.seq
        self.seq = []
        sub = "N" if delete==False else []
        # self.mask is provided as listified string in il_intra
        if len(self.mask) > 0:
            for mask in self.mask:
                maskseq = list((mask[1]-mask[0]) * sub)
                bases_masked += len(maskseq)
                self.inserted[mask[0]:mask[1]] = maskseq
        else:
            pass
        if not logger:
            print(
                f"masked {len(self.mask)} sequences with {bases_masked} bases"
                f" in {self.header}"
                )
        else:
            logger.vlprint(
                f"masked {len(self.mask)} sequences with {bases_masked} bases"
                f" in {self.header}", 2
                )
        self.inserted = "".join(self.inserted)
        return len(self.mask)

    def write_out(self, outfile):
        outfile.write(self.header + "\n")
        outfile.write(self.inserted + "\n")


if __name__ == "__main__":
    args = get_arguments()
    accept = args.acceptor
    acc_fname = os.path.split(accept)[1]
    integr = args.integration
    int_fname = os.path.split(integr)[1]
    number = int(args.number)
    targeted = args.targeted
    time = strftime("%Y-%m-%d_%H-%M", localtime())
    out_name = f"{acc_fname}_{number}xints_{int_fname}_{time}"

    start = perf_counter()
    contigs, total_bp, n_regions = parse_acceptor_seq(accept)

    if targeted:
        int_dict = load_target_sites(targeted, contigs, out_name)
    else:
        int_dict = pick_int_sites(number, contigs, total_bp, out_name, n_regions)
    insert_ints(int_dict, accept, integr, out_name)
    end = perf_counter()
    print(f"Integration locator required {end-start} seconds to complete")
