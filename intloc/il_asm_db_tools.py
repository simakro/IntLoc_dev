# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
import sys
import gzip
import subprocess
from hashlib import md5

from il_logging import Ilogger
from il_subprocesses import caplog_sp_error


ilog = Ilogger()
ilog.module = __name__

def gen_blastdb(_genome, args, **kwargs):
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
    ilog.vlprint(f"Generating BLAST-DB for {msg}", 1)
    try:
        if ena:
                subprocess.run(
                            [
                            args.makeblastdb,
                            "-in", _genome,
                            "-dbtype", "nucl"
                            ],
                        check=True, capture_output=True
                        )
        else:
            subprocess.run(
                            [
                            args.makeblastdb,
                            "-in", _genome,
                            "-dbtype", "nucl",
                            "-parse_seqids"
                            ],
                        check=True, capture_output=True
                        )
    except subprocess.CalledProcessError as e:
        err_txt = e.stderr.decode('utf-8')
        if "ENA|" in err_txt:
            gen_blastdb(_genome, args, msg=msg, ena=True)
        else:
            caplog_sp_error(e, f"makeblastdb for {msg} failed due to {err_txt}")


def check_for_refdb(genome, args):
    if os.path.exists(genome + '.nsq'):
        ilog.vlprint('Database for reference genome already exists.', 2)
    elif genome.split(".")[-1] == "gz":
        comp_genome = genome
        genome = ".".join(genome.split(".")[:-1])
        with open(genome, "w") as extracted:
            gio = gzip.open(comp_genome)
            for line in gio:
                extracted.write(line.decode("utf-8"))
        gen_blastdb(genome, args)
    else:
        gen_blastdb(genome, args)
    return genome


def cleanup_blastdb(reads_dir):
    """Removes blast-database files created by makeblastdb (but not the 
    read/reference files of course)"""
    ilog.vlprint("Cleaning up BLASTdb files.", 2)

    bd_ext = [
        "ndb",
        "nal",
        "nsq",
        "nhr",
        "nin",
        "nog",
        "nsd",
        "nsi",
        "nos",
        "not",
        "nto",
        "ntf"
        ]
    for entry in os.scandir(reads_dir):
        head, tail = os.path.split(str(entry.path))
        ext = tail.split(".")[-1]
        if ext in bd_ext and entry.is_file():
            os.remove(str(entry.path))
#################################################

# From il_asm_db_tools.py#
#################################################
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


def get_genome_info_blast_tool(genome_db):
    """Get reference info from BLAST DB"""
    
    ilog.vlprint("Retrieving genome info from BLAST DB", 2)
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
    return entries_dict


def scan_ref_genome_tool(genome_db, blastdb_info, standalone=False):
    """Custom scan of reference genome"""
    ilog.vlprint("Scanning reference genome with scan_ref_genome_tool", 2)
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
        blastdb_info[cont_id]["len"]["scanref"] = entries_len[cont_id]
        blastdb_info[cont_id]["descr"]["scanref"] = " ".join(entries_descr[cont_id])
        blastdb_info[cont_id]["hash"]["scanref"] = entries_hash[cont_id].hexdigest()

    stats = [entries_len[tig] for tig in entries_len]
    total_bp = sum(stats)
    ilog.vlprint("Genome reference statistics:", 10)
    ilog.vlprint(f"No. of contigs {len(stats)}", 10)
    ilog.vlprint(f"longest contig {max(stats)}", 10)
    ilog.vlprint(f"Total base count {total_bp}", 10)
    return blastdb_info