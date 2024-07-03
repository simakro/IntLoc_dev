# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import sys
import gzip
import shutil
import os

def prune_feature_table(full_table):

    if full_table.split(".")[-1] == "gz":
        mode = "compressed"
        print(mode)
        ff_tab = gzip.open(full_table, "rb")
    else:
        mode = "decompressed"
        print(mode)
        ff_tab = open(full_table, "r")


    if mode == "compressed":
        with gzip.open(full_table + ".slim", "wb") as slim:
            for line in ff_tab:
                # line = line.decode("utf-8")
                ls = line.split(b"\t")
                ignore = [b"mRNA", b"CDS", b"misc_RNA", b"ncRNA", b"tRNA"]
                if ls[0] not in ignore:
                    slim.write(line)
    else:
        with open(full_table + ".slim", "w") as slim:
            for line in ff_tab:
                line = line.replace(",", ";")
                ls = line.split("\t")
                ignore = ["mRNA", "CDS", "misc_RNA", "ncRNA", "tRNA"]
                if ls[0] not in ignore:
                    slim.write(line)
    
    return mode


def compress_and_rename(full_table):
        slim_ver = full_table + ".slim"
        rootname = ".".join(slim_ver.split(".")[:-2])
        newname = rootname + "_slim.txt"
        os.rename(slim_ver, newname)
        with open(newname, "rb") as fin:
            with gzip.open(newname + ".gz", "wb") as fout:
                shutil.copyfileobj(fin, fout)


if __name__ == "__main__":
    full_table = sys.argv[1]
    mode = prune_feature_table(full_table)
    if mode == "decompressed":
        compress_and_rename(full_table)
