# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
import sys
import json
import gzip
import tarfile
import ftplib
import hashlib

from il_logging import Ilogger

ilog = Ilogger()
ilog.module = __name__

def chk_md5_sum(file):
    md5_sum = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_sum.update(chunk)
    return md5_sum.hexdigest()


def get_asm_dwnl_params(species, exact_match=False):
    pkg_dir = os.path.split(os.path.abspath(__file__))[0]
    dwnl_db = os.path.join(pkg_dir, "il_resources", "genome_dwnl_db.json")
    with open(dwnl_db, "r") as db:
        data = db.read()
        ftp_dict = json.loads(data)
    ftp_path = ftp_dict[species]
    ftp_ncbi = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp_ncbi.login()
    ftp_ncbi.cwd(f"{ftp_path}")
    ftp_dir = ftp_ncbi.mlsd()
    latest_asm = []
    for entry in ftp_dir:
        if entry[0].startswith("GCF_"):
            latest_asm.append(entry[0])
    if len(latest_asm) > 1:
        ilog.vlprint(
            f"INFO: More than one assembly file was found in latest assembly"
            f" folder: {latest_asm}. Using {latest_asm[0]}", 1
            )
    if not exact_match:
        genome = latest_asm[0]
    else:
        if exact_match in latest_asm:
            genome = exact_match
        else:
            ilog.vlprint(
                f"Assembly {exact_match} was not found in {ftp_path}. Check "
                f"spelling or use only {species} to get latest assembly.", 0
                )
            sys.exit()

    return ftp_path, genome


def get_genome_asm(
    ftp_path,
    genome,
    pkg_dir=os.path.split(__file__)[0],
    res_dir="il_resources",
    asm_dir="genome_assemblies"
    ):
    gen_file_name = f"{genome}_genomic.fna.gz"
    feat_tbl_file = f"{genome}_feature_table.txt.gz"
    gen_path = os.path.join(pkg_dir, res_dir, asm_dir, genome)
    try:
        os.mkdir(gen_path)
    except FileExistsError as fe:
        ilog.vlprint(str(fe), 1)
        ilog.vlprint(
            f"WARNING: Directory for assembly {gen_path} already exists. "
            f"Intloc will check if the downloaded files in this directory "
            f"appear to be in order and attempt to run integration search with "
            f"them. If this fails, please delete the assembly folder and "
            f"download it again."
            , 0
            )

    gen_fa = os.path.join(gen_path, gen_file_name)
    ftl_txt = os.path.join(gen_path, feat_tbl_file)
    mdt = os.path.join(gen_path, "md5checksums.txt")

    ftp_ncbi = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp_ncbi.login()
    ftp_ncbi.cwd(f"{ftp_path}/{genome}")
    with open(gen_fa, "wb") as gem, open(ftl_txt, "wb") as ftl, open(mdt, "wb") as chk:
        ftp_ncbi.retrbinary(f"RETR {genome}_genomic.fna.gz", gem.write)
        ftp_ncbi.retrbinary(f"RETR {genome}_feature_table.txt.gz", ftl.write)
        ftp_ncbi.retrbinary("RETR md5checksums.txt", chk.write)

    with open(mdt, "r") as mdf:
        downloads = {
            f"{genome}_genomic.fna.gz": "",
            f"{genome}_feature_table.txt.gz": ""
            }
        for line in mdf:
            ls = line.strip().split("  ")
            file = ls[1].split("/")[1]
            if file in downloads:
                md5sum = ls[0]
                downloads[file] = md5sum
    dwnl_files = os.scandir(gen_path)

    dwnl_success = []
    for entry in dwnl_files:
        if entry.name in downloads:
            entr_md5 = chk_md5_sum(entry.path)
            if entr_md5 == downloads[entry.name]:
                ilog.vlprint(f"{entry} was downloaded correctly.", 0)
                dwnl_success.append(True)
            else:
                ilog.vlprint(
                    f"{entry} md5sum did not match. Downloaded file may be"
                    f" corrupted.", 0
                    )
                dwnl_success.append(False)

    ftp_ncbi.quit()
    return all(dwnl_success), gen_path


def extract_gz(gzipped_file):
    unzipped_file = ".".join(gzipped_file.split(".")[:-1])
    with open(unzipped_file, "w") as extracted:
                    gio = gzip.open(gzipped_file)
                    for line in gio:
                        extracted.write(line.decode("utf-8"))
    os.remove(gzipped_file)
    return unzipped_file


def run_getasm(species):
    if "=" in species:
        species, asm_name = species.split("=")
    else:
        asm_name = False
    ftp_path, genome = get_asm_dwnl_params(species, exact_match=asm_name)
    dwnl_success, genome_dir = get_genome_asm(ftp_path, genome)
    genome_path = os.path.join(genome_dir, f"{genome}_genomic.fna.gz")
    genome_path = extract_gz(genome_path)
    # genome_gz = genome_path
    # genome_path = ".".join(genome_gz.split(".")[:-1])
    # with open(genome_path, "w") as extracted:
    #                 gio = gzip.open(genome_gz)
    #                 for line in gio:
    #                     extracted.write(line.decode("utf-8"))
    # os.remove(genome_gz)
    return dwnl_success, genome_path


def get_win_blast(
    pkg_dir=os.path.split(__file__)[0],
    deps_dir="il_deps"
    ):
    deps_path = os.path.join(pkg_dir, deps_dir)
    blast_path = os.path.join(deps_path, "ncbi-blast")
    # from pathlib import Path
    # blastn_bin = list(Path(deps_path).glob(os.path.join("**", "blastn.exe")))
    # mkbldb_bin = list(Path(deps_path).glob(os.path.join("**", "makeblastdb.exe")))


    # try:
    #     ilog.vlprint("Checking if blast dependency dir exists", 0)
    #     # ilog.vlprint("Checking if autodownloaded blast dependency is available", 0)
    #     blast_winbin = os.path.join(blast_path, "*", "bin", "blastn.exe")
    #     avail_bin = [f.path for f in os.scandir(blast_path)] 
    #     scan = [f.path for f in os.scandir(blast_path)]
    #     # print(scan)
    # except:
    #     ilog.vlprint("Autodownloaded blast dependency is not available", 0)



    try:
        ilog.vlprint("Generating ncbi-blast folder in il_deps dir", 0)
        os.makedirs(blast_path, exist_ok=True)
    except FileExistsError as fe:
        ilog.vlprint(str(fe), 1)
        ilog.vlprint(
            f"WARNING: Directory for dependency blast already exists. "
            f"Intloc will check if the downloaded files in this directory "
            f"appear to be in order and attempt to run integration search with "
            f"them. If this fails, please delete the folder {blast_path} and "
            f"download it again."
            , 0
            )
        

    # gen_fa = os.path.join(gen_path, gen_file_name)
    # ftl_txt = os.path.join(gen_path, feat_tbl_file)
    # mdt = os.path.join(gen_path, "md5checksums.txt")

    ilog.vlprint("Logging into ncbi ftp", 0)
    ftp_ncbi = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp_ncbi.login()
    ftp_path = "blast/executables/blast+/LATEST"
    ilog.vlprint("Navigating to ftp download path", 0)
    ftp_ncbi.cwd(f"{ftp_path}")
    #ftp_ncbi.dir()
    ilog.vlprint("Listing  files available for download", 0)
    file_lst = ftp_ncbi.nlst()
    win_archs = [f for f in file_lst if "x64-win64.tar.gz" in f]
    ilog.vlprint(f"Win-files available for download {win_archs}", 5)
    print(win_archs)
    win_archs = {f.split(".")[-1]: f for f in win_archs}
    blast_tar = os.path.join(blast_path, win_archs["gz"])
    blast_md5 = os.path.join(blast_path, win_archs["md5"])
    print(blast_tar)
    print(blast_md5)

    with open(blast_tar, "wb") as bgz, open(blast_md5, "wb") as chk:
        ftp_ncbi.retrbinary(f"RETR {win_archs['gz']}", bgz.write)
        ftp_ncbi.retrbinary(f"RETR {win_archs['md5']}", chk.write)
    
    print(blast_tar)
    print(blast_md5)
    # with tarfile.open(blast_tar) as tar, tarfile.open(blast_md5) as bmd5:
    #     tar.extractall(blast_path)
    #     bmd5.extractall(blast_path)

    
    tar = tarfile.open(blast_tar)
    tar.extractall(blast_path)
    # bmd5 = tarfile.open(blast_md5)
    # bmd5.extractall(blast_path)
    

    # with open(blast_md5, "r") as mdf:
    #     downloads = {
    #         f"{genome}_genomic.fna.gz": "",
    #         f"{genome}_feature_table.txt.gz": ""
    #         }
    #     for line in mdf:
    #         ls = line.strip().split("  ")
    #         file = ls[1].split("/")[1]
    #         if file in downloads:
    #             md5sum = ls[0]
    #             downloads[file] = md5sum
    # dwnl_files = os.scandir(gen_path)

    # dwnl_success = []
    # for entry in dwnl_files:
    #     if entry.name in downloads:
    #         entr_md5 = chk_md5_sum(entry.path)
    #         if entr_md5 == downloads[entry.name]:
    #             ilog.vlprint(f"{entry} was downloaded correctly.", 0)
    #             dwnl_success.append(True)
    #         else:
    #             ilog.vlprint(
    #                 f"{entry} md5sum did not match. Downloaded file may be"
    #                 f" corrupted.", 0
    #                 )
    #             dwnl_success.append(False)

    ftp_ncbi.quit()
    # return all(dwnl_success), gen_path

if __name__ == "__main__":
    # https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/
    # genome = sys.argv[1]
    # genome = "GCF_000001405.38_GRCh38.p12"

    # import sys
    # species = sys.argv[1]
    # run_getasm(species)

    get_blast_dep()
