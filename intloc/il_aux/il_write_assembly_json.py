# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
import gzip
import sys
import json
import subprocess
from collections import defaultdict
from hashlib import md5

# intloc_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# sys.path.insert(0, intloc_dir)

# from il_logging import Ilogger
# # from integration_locator import gen_blastdb
# from il_chk_conv import check_for_refdb


# ilog = Ilogger()
# ilog.module = __name__

def check_for_refdb(genome):
    if os.path.exists(genome + '.nsq'):
        # ilog.vlprint('Database for reference genome already exists.', 1)
        print('Database for reference genome already exists.')
    elif genome.split(".")[-1] == "gz":
        comp_genome = genome
        genome = ".".join(genome.split(".")[:-1])
        with open(genome, "w") as extracted:
            gio = gzip.open(comp_genome)
            for line in gio:
                extracted.write(line.decode("utf-8"))
        gen_blastdb(genome)
    else:
        gen_blastdb(genome)
    return genome


def gen_blastdb(_genome, **kwargs):
    ena = False
    if "msg" not in kwargs:
        msg="reference genome"
    elif "ena" in kwargs:
        ena = True
        msg="ref genome without -parse_seqids due to ENA seq-id incompatibility"
    else:
        msg=kwargs["msg"]
    workpath, genomepath = os.path.split(_genome)
    os.chmod(workpath, 0o777)
    os.chdir(workpath)
    # ilog.vlprint(f"Generating BLAST-DB for {msg}", 0)
    print(f"Generating BLAST-DB for {msg}", 0)
    try:
        if ena:
                subprocess.run(
                            [
                            "makeblastdb",
                            "-in", _genome,
                            "-dbtype", "nucl"
                            ],
                        check=True, capture_output=True
                        )
        else:
            subprocess.run(
                            [
                            "makeblastdb",
                            "-in", _genome,
                            "-dbtype", "nucl",
                            "-parse_seqids"
                            ],
                        check=True, capture_output=True
                        )
    except subprocess.CalledProcessError as e:
        err_txt = e.stderr.decode('utf-8')
        if "ENA|" in err_txt:
            gen_blastdb(_genome, msg=msg, ena=True)
        else:
            print(f"makeblastdb in il_write_genome_info failed {e}")


# def gen_blastdb(_genome, **kwargs):
#     if "msg" not in kwargs:
#         msg="reference genome"
#     else:
#         msg=kwargs["msg"]
#     workpath, genomepath = os.path.split(_genome)
#     os.chmod(workpath, 0o777)
#     os.chdir(workpath)
#     print(f"Generating BLAST-DB for {msg}", 0)
#     try:
#         subprocess.run(
#             ["makeblastdb", "-in", _genome, "-parse_seqids", "-dbtype", "nucl"],
#             check=True, capture_output=True
#             )
#     except subprocess.CalledProcessError as e:
#         print(f"makeblastdb in il_write_genome_info failed {e}")


def process_tig_id(cont_id):
    # gnl = general db identifier; lcl= local seq identifier
    # ref= NCBI reference sequence; bbs = GenInfo Backbone Id
    # for GenBank, EMBL Data Library and DDBJ:
    # gi = gi|gi-number|gb or emb or dbj|accession|locus
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
    return cont_id


def get_genome_info_blast(genome_db):
    """Get reference info from BLAST DB"""
    
    print("Retrieving genome info from BLAST DB")
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
        espl = entry.strip().split("::$&~::")
        if entry.startswith("gnl"):
            cont_id, descr = espl[1], espl[2]
        else:
            cont_id, descr = espl[0], espl[1]
        cont_len, cont_hash = espl[-2], espl[-1]
        if cont_id.startswith("ENA"):
            ci_spl = cont_id.split(" ")
            cont_id, descr = ci_spl[0], " ".join(ci_spl[1:])
        cont_id = process_tig_id(cont_id)
        entries_dict[cont_id] = {
                                "descr": {
                                    "blastdbcmd": descr,
                                    "scanref": None
                                },
                                "len": {
                                    "blastdbcmd": int(cont_len),
                                    "scanref": None
                                },
                                "hash": {
                                    "blastdbcmd": cont_hash,
                                    "scanref": None
                                }
                            }
    # print("blast entries_dict", entries_dict)
    # import pickle
    # pickle(entries_dict)
    return entries_dict


def scan_ref_genome_dct(genome_db, blastdb_info, standalone=False): # int_dict
    """Custom scan of reference genome"""
    print(
        "scanning reference genome"
        )
    entries_len = {}
    entries_descr = {}
    entries_hash = {}
    curr_cont = ""
    with open(genome_db, "r") as gen:
        for line in gen:
            if line.startswith(">"):
                if len(curr_cont) > 0:
                    entries_hash[curr_cont] = "".join(entries_hash[curr_cont])
                    entries_hash[curr_cont] = md5(entries_hash[curr_cont].encode())
                split_tig = line.strip().replace(">", "").split(" ")
                if line.startswith("gnl"):
                    cont_id =  split_tig[1]
                    descr = split_tig[2:]
                else:
                    cont_id = split_tig[0]
                    descr = split_tig[1:]
                cont_id = process_tig_id(cont_id)
                curr_cont = cont_id
                entries_len[curr_cont] = 0
                entries_hash[curr_cont] = []
                entries_descr[curr_cont] = descr
            else:
                entries_len[curr_cont] += len(line.strip())
                entries_hash[curr_cont].append(line.strip())
        entries_hash[curr_cont] = "".join(entries_hash[curr_cont])
        entries_hash[curr_cont] = entries_hash[curr_cont].upper()
        entries_hash[curr_cont] = md5(entries_hash[curr_cont].encode())
        # entries_hash = {k:"".join(v) for k,v in entries_hash.items()}
        # entries_hash = {k:md5(v.encode()) for k,v in entries_hash.items()}

    if standalone:
        blank_entry = {
                        "descr": {
                            "blastdbcmd": None,
                            "scanref": None
                        },
                        "len": {
                            "blastdbcmd": None,
                            "scanref": None
                        },
                        "hash": {
                            "blastdbcmd": None,
                            "scanref": None
                        }
                    }
        blastdb_info = {tig: blank_entry for tig in entries_len}

    for cont_id in entries_len:
        # split_tig = contig.replace(">", "").split(" ")
        # if contig.startswith("gnl"):
        #     cont_id, descr =  split_tig[1], split_tig[2]
        # else:
        #     cont_id, descr = split_tig[0], split_tig[1]
        # descr = " ".join(descr)
        blastdb_info[cont_id]["len"]["scanref"] = entries_len[cont_id]
        blastdb_info[cont_id]["descr"]["scanref"] = " ".join(entries_descr[cont_id])
        blastdb_info[cont_id]["hash"]["scanref"] = entries_hash[cont_id].hexdigest()

    # print("blastdb + scanref info", blastdb_info)

    stats = [entries_len[tig] for tig in entries_len]
    total_bp = sum(stats)
    # ilog.vlprint("Genome reference statistics:", 10, outdir=outdir)
    # ilog.vlprint(f"No. of contigs {len(stats)}", 10, outdir=outdir)
    # ilog.vlprint(f"longest contig {max(stats)}", 10, outdir=outdir)
    # ilog.vlprint(f"Total base count {total_bp}", 10, outdir=outdir)
    print("Genome reference statistics:")
    print(f"No. of contigs {len(stats)}")
    print(f"longest contig {max(stats)}")
    print(f"Total base count {total_bp}")
    return blastdb_info


def write_json(entries_dict, rg_id, species, gen_res_dbfile):
    with open(gen_res_dbfile, "r") as grdb:
        data = grdb.read()
        jdat = json.loads(data)
        jdat = defaultdict(dict, jdat)
    print("jdat", jdat)
    # could add ["source"] or ["provider"] (e.g. NCBI, ensembl etc.) after 
    # ["rg_id"] to account for potentially different tig names etc. 
    jdat[species][rg_id] = entries_dict 
    # new_entry = json.load(entries_dict)
    print("new_entry added", jdat)

    # out_new = gen_res_dbfile + "_new"
    # with open(out_new, "w") as outf:
    
    for species in jdat:
        jdat[species] =  dict(sorted(jdat[species].items(), key=lambda x: x[0]))
    jdat = dict(sorted(jdat.items(), key=lambda x: x[0]))
    with open(gen_res_dbfile, "w") as outf:
        out_data = json.dumps(jdat, indent=2)
        outf.write(out_data)

def compare_assembly_ident(species, q_assembly, genome_ident_db, data_source="blastdbcmd"):
    with open(genome_ident_db, "r") as db:
        data = db.read()
        jdat = json.loads(data)

    for asm in jdat[species]:
        if asm != q_assembly:
            for entry in jdat[species][q_assembly]:
                # if entry.startswith("NC_"):
                try:
                    comp_len = "" if jdat[species][asm][entry]["len"][data_source] == jdat[species][q_assembly][entry]["len"][data_source] else "not "
                    print(f"{asm} {entry} length is {comp_len}identical to {q_assembly}")
                    comp_len = "" if jdat[species][asm][entry]["hash"][data_source] == jdat[species][q_assembly][entry]["hash"][data_source] else "not "
                    print(f"{asm} {entry} hash value is {comp_len}identical to {q_assembly}")
                except KeyError:
                    print(f"Contig {entry} is not present in {asm}")





if __name__ == "__main__":
    # path to reference genome
    ref_gen = sys.argv[1]
    # reference genome id (key for assembly in json/dict)
    asm_id = sys.argv[2]
    # the genome assembly identification database json file
    gen_res_dbfile = sys.argv[3]
    # formal species name e.g. Homo sapiens, Mus musculus etc. 
    species = sys.argv[4]

    check_for_refdb(ref_gen)
    blast_entries = get_genome_info_blast(ref_gen)
    entries_dict = scan_ref_genome_dct(ref_gen, blast_entries)
    write_json(entries_dict, asm_id, species, gen_res_dbfile)

    # print(entries_dict)
    

    # # for comparison of individual assmeblies in db
    # species = sys.argv[1]
    # assembly = sys.argv[2]
    # genome_ident_db = sys.argv[3]
    # compare_assembly_ident(species, assembly, genome_ident_db)
