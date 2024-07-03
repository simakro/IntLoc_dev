# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
import sys

intloc_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, intloc_dir)

from il_get import get_genome_asm

def download_asm_lst(asm_lst, out_dir, res_dir, asm_dir):
    """Download a list of assemlies from NCBI"""
    for genome in asm_lst:
        get_genome_asm(
            genome,
            pkg_dir=out_dir,
            res_dir=res_dir,
            asm_dir=asm_dir
            )


if __name__ == "__main__":

    asm_lst = [
        # "GCF_000001405.30_GRCh38.p4",
        # "GCF_000001405.31_GRCh38.p5",
        # "GCF_000001405.32_GRCh38.p6",
        # "GCF_000001405.33_GRCh38.p7",
        # "GCF_000001405.34_GRCh38.p8",
        # "GCF_000001405.35_GRCh38.p9",
        # "GCF_000001405.36_GRCh38.p10",
        "GCF_000001405.37_GRCh38.p11",
        ]
    # out_dir = "/Users/user/Documents/reference_genomes"
    out_dir = sys.argv[1]
    res_dir = "human"
    asm_dir = "script_downloads"

    download_asm_lst(asm_lst, out_dir, res_dir, asm_dir)

