# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

def read_integrator_report(report):
    with open(report, "r") as rep:
        int_sites = dict()
        csv_line = False
        int_len = int()
        for line in rep:
            if line.startswith("Integrating sequence length"):
                int_len = int(line.strip().split(": ")[1])
            elif len(line.strip()) > 0:
                try:
                    int(line[0])
                    csv_line = True
                except ValueError as ve:
                    csv_line = False

                if csv_line:
                    int_num, chrom, int_site = line.strip().split(",")
                    int_num, chrom, int_site = int(int_num), chrom.strip(), int(int_site.strip())
                    if chrom in int_sites:
                        int_sites[chrom].append((int_num, int_site))
                    else: 
                        int_sites[chrom] = [(int_num, int_site)]
    print(int_sites, int_len)
    return int_sites, int_len


def crop_sequence(int_sites, int_len, margins, fa_file):
    fasta = open(fa_file, "r")
    with open("cropped_sequences.fasta", "w") as out:
        cur_chrom = ""
        cur_seq = ""
        ct = 0
        for line in fasta:
            if line.startswith(">"):
                ct += 1
                if ct == 1:
                    cur_chrom = line.strip()
                else:
                    if cur_chrom in int_sites:
                        chrom_name = cur_chrom.split(">")[1]
                        for int_num, int_site in int_sites[cur_chrom]:
                            out.write(f">Integration_{int_num}_in_{chrom_name}_crop_with_{margins}_margins\n")
                            out.write(cur_seq[int_site-1-margins: int_site-1+int_len-1+margins] + "\n")
                    else:
                        pass
                    cur_chrom = line.strip()
                    cur_seq = ""
            else:
                if cur_chrom in int_sites:
                    cur_seq = cur_seq + line.strip()
                else:
                    pass


if __name__ == "__main__":
    report = "./integrator_report.csv"
    margins = 2500
    fasta = "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_ints2.fna"
    int_sites, int_len = read_integrator_report(report)
    crop_sequence(int_sites, int_len, margins, fasta)
