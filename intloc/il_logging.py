# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import sys
import json
import os
from time import strftime, localtime


class Ilogger():
    
    def __init__(self):
        self.outdir = str()
        self.verbosity = int()
        self.module = str()


    def report_module(self):
        print(self.module) 

    def vlprint(self, report, level, logging=True, alt_log=False):
        """Verbosity level print function. For debugging, logging and output."""
        if not logging:
            target = sys.stdout
        elif alt_log:
            target = open(os.path.join(self.outdir, alt_log), "a")
        else:
            target = open(os.path.join(self.outdir, "intlog.log"), "a")

        try:
            if logging:
                time_stamp = strftime("%Y-%m-%d %H:%M:%S", localtime())
                print(f"{time_stamp} {self.module} {report}", file=target, flush=True)
                if level <= int(self.verbosity):
                    print(report, file=sys.stdout, flush=True)
            elif level <= int(self.verbosity):
                print(report, file=target, flush=True)
            else:
                pass
        except NameError:
            pass

        if logging:
            target.close()


def initialize_loggers(args):
    """Set correct verbosity level for loggers of all intloc modules"""
    with open(os.path.join(args.res_dir, "log_init_modules.json"), "r") as data:
        init_logger = json.load(data)
        for mod in init_logger["modules"]:
            sys.modules[mod].ilog.verbosity = args.verbosity