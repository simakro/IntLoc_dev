import os
import sys
from time import perf_counter
from il_logging import Ilogger

ilog = Ilogger()
ilog.module = __name__

def retrieve_reads(filtered_blast_report, reads, read_dir=False, prefix="Integration"):
    """Retrieve candidate reads containing integration from sequencing data"""
    ilog.vlprint('Retrieve candidate reads', 1, logging=True)
    if os.path.exists('Integration_candidate_reads.fasta') and not read_dir:
        ilog.vlprint('Integration_candidate_reads.fasta already exists. Skipping retrieval of reads from source fasta', 1,
                logging=True)
    else:
        search_list = []
        with open(filtered_blast_report, 'r') as fbr:
            for line in fbr:
                if line.startswith('RNAME'):
                    pass
                else:
                    l = line.split(',')
                    header = ">" + l[2]
                    if header in search_list:
                        pass
                    else:
                        search_list.append(header)
    
        print("Length search list", len(search_list))
        print(search_list)
        cand_read_path = os.path.join(os.getcwd(), f"{prefix}_candidate_reads.fasta")
        match = False
        int_reads = open(cand_read_path, "a")
        with open(reads, 'r') as fileo:
            for line in fileo:
                if match is False:
                    if line.startswith('>'):
                        line_tmp = line.strip()
                        hls = line_tmp.split(" ")[0]
                        if hls in search_list:
                            # line_tmp = line
                            not_accepted = [",", ":", "|", ";", "\\", "/"]
                            lt_scan = [sym for sym in line_tmp if sym in not_accepted]
                            if len(lt_scan) > 0:
                                for sym in lt_scan:
                                    line_tmp = line_tmp.replace(sym, " ")
                            line_tmp = line_tmp.strip().split(" ")
                            # if line_tmp[0] == hls:
                            print(line_tmp[0] + "\n", end='', file=int_reads)
                            match = True
                    else:
                        pass
                elif match:
                    print(line, end='', file=int_reads)
                    match = False
        int_reads.close()
    return cand_read_path


if __name__ == "__main__":
    filtered_blast_report = sys.argv[1]
    reads = sys.argv[2]
    start = perf_counter()
    retrieve_reads(filtered_blast_report, reads)
    end = perf_counter()
    print(f"Required {end-start} seconds")