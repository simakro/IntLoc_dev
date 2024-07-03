# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import gzip
from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

class IntFeatures:

    def __init__(
        self, ID, name, chr, loc, loc_acc, re_supp, mms, non_int_rs, cov, IMR, SSBR, rnames
        ):
        self.ID = ID
        self.name = name
        self.chr = chr
        self.loc = loc
        self.loc_acc = loc_acc
        self.re_supp = re_supp
        self.IMR = IMR
        self.SSBR = SSBR
        self.mms = mms
        self.non_int_rs = non_int_rs
        self.cov = cov
        self.read_names = rnames
        self.features = []


def il_features(features, int_rep_csv, prefix=""):
    ilog.vlprint(
        "Identifying genomic feature information for integration sites", 3
        )
    integrations = []
    carr_chr = {}

    with open(int_rep_csv, "r") as int_rep_csv:
        for line in int_rep_csv:
            if line.startswith("ID"):
                pass
            else:
                ls = line.split(",")
                int_feat = IntFeatures(*ls)
                integrations.append(int_feat)
                if ls[2].strip() in carr_chr:
                    carr_chr[ls[2].strip()].append(int(ls[3].strip()))
                else:
                    carr_chr[ls[2].strip()] = [int(ls[3])]

    for line in features:
        line = line.decode("utf-8")
        ls = line.split("\t")
        if ls[6] in carr_chr:
            for int_loc in carr_chr[ls[6]]:
                if int(ls[7]) <= int_loc:
                    if int_loc <= int(ls[8]):
                        ignore = ["mRNA", "CDS", "misc_RNA", "ncRNA"]
                        if ls[0] not in ignore:
                            for i in integrations:
                                if int(i.loc) == int(int_loc) and i.chr.strip() == ls[6].strip():
                                    i.features.append(line)

    with open(f"{prefix}int_features.csv", "w") as feat:
        feat.write("Integration, genomic features\n")
        for i in integrations:
            i.features = [feature.strip() for feature in i.features]
            feat_str = "::".join(i.features)
            if len(i.features) > 0:
                feat.write(f"{i.name},{feat_str}\n")
            else:
                feat.write(f"{i.name},no feature reported\n")


if __name__ == "__main__":
    import sys

    _features = gzip.open(sys.argv[1], "r")
    _int_rep_csv = sys.argv[2]

    il_features(_features, _int_rep_csv)
