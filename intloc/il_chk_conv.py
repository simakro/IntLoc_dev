# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.


import argparse
import os
import sys
import json
import math
# import gzip
import subprocess
from time import localtime, strftime, perf_counter  #, sleep
from collections import defaultdict, Counter

# from il_aux.il_write_assembly_json import check_for_refdb #, gen_blastdb
# from il_aux.il_write_assembly_json import get_genome_info_blast, scan_ref_genome
from il_asm_db_tools import check_for_refdb #, gen_blastdb
from il_asm_db_tools import get_genome_info_blast_tool, scan_ref_genome_tool
from il_get import get_win_blast
from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

# def initialize_chkconv_logger(args):
#     ilog.verbosity = args.verbosity

def report_time(module, end=False, alt_log=False):
    """Get time run was started"""
    time_str = strftime("%Y-%m-%d %H:%M:%S", localtime())
    perf_count = perf_counter()
    if not end:
        ilog.vlprint(
            f'\033[1mStarting {module} at: {time_str} \033[0m', 
            0, logging=False
            )
        ilog.vlprint(f'Starting {module} at: {time_str}', 10, alt_log=alt_log)
    else:
        ilog.vlprint(
        f'\033[1m{module} finished at: {time_str}\033[0m',
        0, logging=False
        )
        ilog.vlprint(f'{module} finished at: {time_str}', 10, alt_log=alt_log)
    return time_str, perf_count


def determine_num_threads(args):
    """Determine number of threads to use and adjust args.cores accordingly"""
    avail_cpu = os.cpu_count()
    if args.cores:
        thr = args.cores
        thr = thr if int(thr)<=avail_cpu else str(avail_cpu)
    else:
        thr = avail_cpu - 2
        thr = str(thr) if thr > 1 else "2"
    args.cores = thr
    return thr


def get_full_species_name(args):
    with open(os.path.join(args.res_dir, "spec_name_conv.json"), "r") as sncf:
        data = sncf.read()
        spec_name_conv = json.loads(data)
    try:
        species = spec_name_conv[args.species]
    except KeyError as ke:
        ilog.vlprint(f"WARNING: KeyError {ke}. Wrong species synonym.", 0)
        ilog.vlprint(
            "WARNING: Species specific ideogram and annotations could not"
            " be generated, due to invalid species parameter. Please, check"
            " spelling. E.g. for human genome use either \"human\" or"
            " \"Hs\" or \"Homo sapiens\"", 0
            )
    return species


# def ranged_typechk(min, max, in_type, out_type=False):

#     def ranged_num_type(arg):
#         """Type and range checker for numeric (except complex) cli arguments"""

#         try:
#             arg = eval(arg)
#         except:
#             raise argparse.ArgumentTypeError(f"value has to be a number")
        
#         accepted = [int, float] # complex doesn't make sense here
#         if not type(arg) in accepted:
#             raise argparse.ArgumentTypeError(f"value has to be numeric")
        
#         if type(arg) != in_type:
#             if in_type==int:
#                 raise argparse.ArgumentTypeError(
#                     f"value has to be an integer"
#                     )
#             elif in_type==float:
#                 arg = in_type(arg)
#             else:
#                 raise argparse.ArgumentTypeError(
#                     f"in_type arg must be int or float"
#                     )

#         if not min <= arg <= max:
#             raise argparse.ArgumentTypeError(
#                 f"value has to be within the range of {min} and {max}"
#                 )
        
#         if out_type:
#             arg = out_type(arg)

#         return arg

#     return ranged_num_type


def ranged_typechk(min, max, in_type, out_type=False, alt=False):

    def ranged_num_type(arg):
        """Type and range checker for numeric (except complex) cli arguments"""

        try:
            arg = eval(arg)
        except:
            raise argparse.ArgumentTypeError(f"value has to be a number")
        
        accepted = [int, float] # complex doesn't make sense here
        if alt:
            accepted.append(alt)
        if not type(arg) in accepted:
            raise argparse.ArgumentTypeError(f"value has to be numeric")
        
        if alt==type(arg):
            pass
        else:
            if type(arg) != in_type:
                if in_type==int:
                    raise argparse.ArgumentTypeError(
                        f"value has to be an integer"
                        )
                elif in_type==float:
                    arg = in_type(arg)
                else:
                    raise argparse.ArgumentTypeError(
                        f"in_type arg must be int or float"
                        )

            if not min <= arg <= max:
                raise argparse.ArgumentTypeError(
                    f"value has to be within the range of {min} and {max}"
                    )
            
            if out_type:
                arg = out_type(arg)

        return arg

    return ranged_num_type


def boolmix_typechk(_type):

    def boolstr_type(arg):
        "Type and range checker for numeric cli arguments"
        arg = eval(arg)
        if arg:
            try:
                arg = _type(arg)
            except ValueError:
                raise argparse.ArgumentTypeError(
                    f"Argument has to be of type {_type} or False"
                    )
        return arg

    return boolstr_type


def informed_round(val, info):
    if info > val:
        return math.ceil(val)
    elif info < val:
        return math.floor(val)
    else:
        return round(val)


def required_dependencies(args):
    req_deps = []
    fac_deps = []
    if any([args.intra, args.quick_mm2]):
        req_deps.append("minimap2")
        fac_deps.extend(["paftools.js", "samtools"])
    else:
        # standard and polyclonal mode
        req_deps.extend(["blastn", "makeblastdb", "minimap2"])
        if args.skip_eval:
            req_deps.pop("minimap2")
    return req_deps, fac_deps


def test_samtools_functionality(args):
    """More thorough test of samtools funcionality"""
    ilog.vlprint("Performing samtools fasta dump functionality test", 4)
    ilog.vlprint("Testing if samtools is working", 2)
    test_sam = os.path.join(args.pkg_dir, "il_resources", "test_sam.sam")
    test_out = os.path.join(args.pkg_dir, "il_resources", "test_sam.out")
    sam_test = subprocess.run(
                            [
                                "samtools",
                                "fasta",
                                test_sam,
                                "-F", "4,256,512,2048"
                            ],
                            # stdout=subprocess.PIPE
                            check=True, capture_output=True
                            )
    with open(test_out, "w") as out:
        out.write(sam_test.stdout.decode())
    with open(test_out, "r") as out:
        data = out.read()
    ilog.vlprint(
        f"sam fasta dump contains {len(data)} characters. Expected 43944", 4
        )
    os.remove(test_out)


def chk_minimap_version(mm_ver, args):
    avail = mm_ver.split("-")[0].split(".")
    # reqmaj, reqmin = args.defaults["req_version minimap2"]
    reqmaj, reqmin = [2,16] # -o arg is available from version 2.16 on
    if type(args.dep_ver) != dict:
        args.dep_ver = {} # required for gui, which transforms to str
    try:
        avail = [eval(v) for v in avail]
        if all([avail[0]>=reqmaj, avail[1]>=reqmin]):
            ilog.vlprint("minimap2 version is sufficiently recent", 3)
            args.dep_ver["minimap2"] = "req met"
        else:
            ilog.vlprint(
                "minimap2 version is outdated and may cause problems"
                , 0
                )
            args.dep_ver["minimap2"] = "outdated"
    except:
        ilog.vlprint(
            "Comparison of required to available minimap2 version failed."
            , 2
            )
        args.dep_ver["minimap2"] = "n.a."

def chk_nonpy_dependency(dep, args):
    """Check for availability and basic functionality of dependency"""
    ilog.vlprint(
        f"Checking availability and basic functionality of {dep} dependency", 4
        )
    dep_unavail_msg = f"Intloc will attempt to finalize this run without {dep}. "
    if dep=="minimap2":
        dep_unavail_msg = dep_unavail_msg + f"IntLoc will still work in standard"\
            f" mode, but following modules/modes will not be available:"\
            f" il_intra (--intra), il_polyclonal (--polyclonal), il_quick_mm2"\
            f" (--quick) and il_evaluate."
    if dep=="blastn" or dep=="makeblastdb":
        dep_unavail_msg = dep_unavail_msg + f"IntLoc will not work in standard"\
            f" mode, but following modes should still work: il_quick_mm2,"\
            f" il_intra."
    if dep=="samtools":
        print("Testing samtools is installed")
        dep_unavail_msg = dep_unavail_msg + f"{dep} is a facultative dependenc"\
            f"y in quick and intra modes. However, intloc can substitute"\
            f" for {dep} with internal functions. This may affect runtime"\
            f" of these two modes but will yield equivalent results." 
    if dep =="paftools.js":
        dep_unavail_msg = dep_unavail_msg + f"{dep} is a facultative depend"\
            f"ency in quick and intra modes. However, intloc can substitute "\
            f" for {dep} with internal functions. This may affect runtime"\
            f" of these two modes but will yield equivalent results." 
    test_cmds = {
        "minimap2": "--version",
        "blastn": "-version",
        "makeblastdb": "-version",
        "samtools": "--help",
        "paftools.js": "version"
        }
    cmd = test_cmds[dep]
    try:
        dep_ver = subprocess.run([dep, cmd], stdout=subprocess.PIPE)
        dep_ver = dep_ver.stdout.replace(
                                        b',', b';'
                                    ).replace(
                                            b'\n', b''
                                        ).replace(
                                                b'\r', b''
                                            ).decode("utf-8")
        ilog.vlprint(f"{dep} installed and responsive", 0)
        ilog.vlprint(dep_ver, 2)
        args.__dict__[dep.split(".")[0]] = dep
        if dep=="minimap2":
            chk_minimap_version(dep_ver, args)
    except subprocess.CalledProcessError as ce:
        ilog.vlprint(f"{dep} check CalledProcessError: {ce}", 1)
        ilog.vlprint(
            f"The {dep} executable is installed but appears not"
            f" to be working properly. Check logs for more info.", 1
            )
        print(
            f"The {dep} executable is installed but appears not"
            f" to be working properly. Check logs for more info.", 1
            )
        ilog.vlprint(dep_unavail_msg, 1)
        args.__dict__[dep.split(".")[0]] = False
    except FileNotFoundError as fe:
        ilog.vlprint(f"Dependency check {dep}: FileNotFoundError: {fe}", 1)
        pf = sys.platform
        if all([pf=="win32", dep=="minimap2"]):
            try:
                ilog.vlprint(
                    f"Attempting to use precompiled minimap2 bin"
                    f" packaged into integration_locator", 1
                    )
                il_mm2_bin = os.path.join(
                    args.pkg_dir,
                    "il_deps",
                    "minimap2_v2.7-r1122_win10_64",
                    "minimap2"
                    )
                subprocess.run([il_mm2_bin, "--version"])
                args.__dict__[dep] = il_mm2_bin
            except:
                ilog.vlprint(
                    f"Usage of precompiled intloc minimap2 bin failed. Try to"
                    f" compile minimap2 from source on your system and add to"
                    f" PATH or contact the intloc maintainer.", 1
                    )
                ilog.vlprint(dep_unavail_msg, 1)
                args.__dict__[dep] = False
        # elif pf=="win32" and dep=="blastn" or dep=="makeblastdb":
        elif all([pf=="win32", dep=="blastn"]) or all([pf=="win32", dep=="makeblastdb"]):
                deps_path = os.path.join(args.pkg_dir, "il_deps")
                # blast_path = os.path.join(deps_path, "ncbi-blast")
                from pathlib import Path
                blastn_bin = list(Path(deps_path).glob(os.path.join("**", "blastn.exe")))
                mkbldb_bin = list(Path(deps_path).glob(os.path.join("**", "makeblastdb.exe")))
                if any(blastn_bin) and any(mkbldb_bin):
                    ilog.vlprint(
                        "Autodownloaded blastn and makeblastdb binaries for Win are available", 3
                        )
                    ilog.vlprint(
                        f"Already downloaded paths {blastn_bin[0]} of type {type(blastn_bin[0])}", 3
                        )
                    ilog.vlprint(
                        f"Already downloaded paths {mkbldb_bin[0]} of type {type(blastn_bin[0])}", 3
                        )
                    args.__dict__["blastn"] = str(blastn_bin[0])
                    args.__dict__["makeblastdb"] = str(mkbldb_bin[0])
                else:
                    ilog.vlprint(
                        "No blastn and makeblastdb binaries for Win available."\
                        " Attempting automatic download from ncbi.", 3
                        )
                    try:
                        get_win_blast()
                        ilog.vlprint(
                            f'Newly downloaded paths '
                            f'{list(Path(deps_path).glob(os.path.join("**", "blastn.exe")))[0]}'
                            f' of type '
                            f'{type(list(Path(deps_path).glob(os.path.join("**", "blastn.exe")))[0])}',
                            3
                            )
                        ilog.vlprint(
                            f'Newly downloaded paths '
                            f'{list(Path(deps_path).glob(os.path.join("**", "makeblastdb.exe")))[0]}'
                            f' of type '
                            f'{type(list(Path(deps_path).glob(os.path.join("**", "makeblastdb.exe")))[0])}',
                            3
                            )
                        args.__dict__["blastn"] = str(
                                list(
                                    Path(deps_path).glob(
                                        os.path.join("**", "blastn.exe")
                                        )
                                    )[0]
                                )
                        args.__dict__["makeblastdb"] = str(
                                list(
                                    Path(deps_path).glob(
                                        os.path.join("**", "makeblastdb.exe")
                                        )
                                    )[0]
                                )
                    except:
                        ilog.vlprint(
                            "Automated download of blast binaries failed", 1
                            )
                        ilog.vlprint(dep_unavail_msg, 1)
        # elif pf=="win32" and dep=="samtools" or dep=="paftools.js":
        elif all([pf=="win32", dep=="samtools"]) or all([pf=="win32", dep=="paftools.js"]):
            ilog.vlprint(dep_unavail_msg, 1)
            ilog.vlprint(f"Running intloc on Windows without {dep}", 1)
            args.__dict__[dep.split(".")[0]] = False
        # elif pf=="darwin" and dep=="samtools" or dep=="paftools.js":
        elif all([pf=="darwin", dep=="samtools"]) or all([pf=="darwin", dep=="paftools.js"]):
            ilog.vlprint(dep_unavail_msg, 1)
            ilog.vlprint(
                f"Running intloc on MacOS without {dep}. This will work, "
                f"and won't affect result quality, but you can easily install"
                f" it via conda or get binaries from github etc. if you like."
                , 1
                )
            args.__dict__[dep.split(".")[0]] = False
        # elif all([pf=="linux", dep=="paftools.js"]):
        #     ilog.vlprint(dep_unavail_msg, 1)
        #     ilog.vlprint(
        #         f"Running intloc on Linux without {dep}. This will work, "
        #         f"and won't affect result quality, but you can easily install"
        #         f" it via conda or get binaries from github etc. if you like."
        #         , 1
        #         )
        #     args.__dict__[dep.split(".")[0]] = False
        else:
            bin_url = {
                "minimap2": [
                    "https://github.com/lh3/minimap2/releases",
                    "Linux and MacOs"
                    ],
                "paftools.js": [
                    "https://github.com/lh3/minimap2/releases",
                    "Linux and MacOs"
                    ],
                "blastn": [
                    "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/",
                    "Windows, Linux and MacOs"
                    ],
                "makeblastdb": [
                    "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/",
                    "Windows, Linux and MacOs"
                    ],
                "samtools": [
                    "https://github.com/samtools/samtools/releases",
                    "Linux and MacOs"
                    ],
                }
            ilog.vlprint(
                        f"The {dep} command failed because the {dep} binary"
                        f" was not found. If you already have a binary, add its"
                        f" location to PATH, or if you installed through conda or"
                        f" similar activate the respective environment. On linux"
                        f" and MacOS you can install {dep} through conda. Precompiled"
                        f" binaries of {dep} for {bin_url[dep][1]} systems can also"
                        f" be downloaded from: {bin_url[dep][0]}", 1
                        )
            ilog.vlprint(dep_unavail_msg, 1)
            args.__dict__[dep] = False
    
    if dep=="samtools":
        if args.samtools:
            try:
                test_samtools_functionality(args)
            except:
                print("samtools is not working set args.samtools to False")
                args.samtools = False
                ilog.vlprint(
                    "Funcionality test of installed samtools failed.", 0
                    )


def chk_compatibility_args(
        incompatible: list,
        module: str,
        args: object,
        opt_msg=[]
        ):
    arg_trans = {
        "intra": "intra",
        "quick": "quick_mm2",
        "polyclonal": "polyclonal",
        "ava_clust": "cluster_reads",
        "reads_dir": "reads_dir",
        "req_int_cont": "ric"
    }
    incompatible_args = [args.__dict__[arg_trans[arg]] for arg in incompatible]
    if any(incompatible_args):
        ddash_args = [f"--{arg}" for arg in incompatible]
        incomp_str = ", ".join(ddash_args[:-1])
        incomp_str = incomp_str + f" and {ddash_args[-1]}"
        ilog.vlprint(
            f"IntLoc {module} is not compatible with {incomp_str} arguments.", 1
            )
        for msg in opt_msg:
            ilog.vlprint(msg, 1)
        sys.exit()


def set_module_defaults(
        module: str,
        parameters: list,
        args: object,
        ):
    arg_trans = {
        "req_int_cont": "ric",
        "min_dist_ints": "min_dist_ints",
    }
    for param in parameters:
        std_def = args.defaults[param]
        aka = arg_trans[param]
        if args.__dict__[aka]==std_def:
            module_def = args.defaults[module][param]
            args.__dict__[aka] = module_def
            ilog.vlprint(
                f"Setting {param} to {module} default {module_def}. If you want"
                f" to adjust this param use any other value (except {std_def})."
                , 0
                )


def get_blast_stats(qfile_path: str, outfile_path: str, ref_genome: str, args):
    ilog.vlprint(f"Retrieving blast statistical parameters", 2)
    subprocess.run([
                args.blastn,
                "-db", ref_genome,
                "-query", qfile_path,
                "-out", outfile_path,
                "-outfmt", "15",
                "-num_threads", "2",
                "-perc_identity", "90",
                "-qcov_hsp_perc", "2"
                ],
                check=True, capture_output=True
                )
    blast_params = {}
    with open(outfile_path, "r") as js:
        jsl = json.load(js)
        ka_stats = jsl["BlastOutput2"][0]["report"]["results"]["search"]["stat"]
        blast_params = jsl["BlastOutput2"][0]["report"]["params"]
        blast_params.update(ka_stats)
        excl = ["eff_space", "hsp_len"]
        blast_params = {k:v for k,v in blast_params.items() if k not in excl}
    return blast_params


def check_input_file(_fasta, check, interact=False, construct_fasta=False, silent=False):
    """Check for format and integrity of fasta file. Quick: probe the first """
    """1000 lines. Extensive: Scan whole file."""
    vb = 1 if not silent else 10
    ilog.vlprint(f"Checking input file: {_fasta} .", vb+1)

    l = 0
    # elab = empty lines at beginning of file
    elab = 0 
    converted = False
    fq_conv = False
    fa_header = False
    new_name = ""
    any_seq = 0
    any_def = 0
    bases = ["a", "t", "g", "c"]
    with open(_fasta, 'r') as fasta:
        corrupt = False
        for line in fasta:
            l += 1
            if any_def == 0 and any_seq == 0:
                if len(line.strip()) == 0:
                    elab += 1
                    ilog.vlprint("Input file starts with empty lines", vb-1)
                    corrupt = True
            if line.startswith(">"):
                any_def += 1
                fa_header = line.strip()
            else:
                if len(line.strip()) > 0:
                    fa_header = False
                if line[0].lower() in bases:
                    any_seq += 1
                    if any_def == 0:
                        ilog.vlprint(
                            "Encountered sequence before first fasta header." 
                            " Input file corrupted.", vb-1
                        )
                        sys.exit()

            if (l - elab) == 1000:
                if check == "quick":
                    break
            elif line.startswith("@"):
                if corrupt:
                    ilog.vlprint(
                    "You provided an irregular fastq file. Intloc can convert"
                    " regular fastq to the preferred fasta input format. Please"
                    " check the input read file/s.", vb-1
                    )
                else:
                    ilog.vlprint(
                        "Incompatible read file input format: Fasta input is "
                        "required, not fastq.", vb-1
                        )
                    if interact:
                        answer = input(
                            "Integration_locator can generate a fasta copy of"
                            " input fastq file. Do you wish to proceed (Y/N)?"
                            ).upper()
                    else:
                        ilog.vlprint(
                            "Integration_locator will generate a fasta copy"
                            " of the input fastq file.", vb
                        )
                        answer = 'Y'
                    corrupt = True
                    if answer == 'N':
                        ilog.vlprint(
                            "Found fastq read input file. User declined"
                            " conversion to fasta. Exiting intloc.", vb-1
                            )
                        sys.exit()
                    elif answer == 'Y':
                        converted = True
                        fq_conv = True
                        try:
                            new_name = fastq_to_fasta(_fasta)
                            corrupt = False
                            ilog.vlprint(
                                "Input fastq read file was converted to fasta.",
                                vb
                                )
                        except:
                            ilog.vlprint(
                                "An unexpected error occurred while attempting"
                                " to generate converted input file. Please"
                                " check input fastq for irregularities", vb-1
                                )
                            raise
                        break
            elif (l - elab) % 2 == 0:
                if line.startswith(">"):
                    ilog.vlprint(
                        "Fasta file is either in Multi-line format (Single-line"
                        " format is required) or contains newlines or other "
                        "entries between sequence records.", vb-1
                        )
                    corrupt = True
                    ilog.vlprint(f"Encountered first error at line {l}", vb-1)
                    break
                else:
                    if len(line.strip()) == 0:
                        corrupt = True
                        if any_def:
                            ilog.vlprint(
                                f"Empty line after defline {fa_header} in"
                                f" line {l}. Missing sequence information or"
                                f" interspersed blank lines. Fasta may be"
                                f" corrupted.", vb-1
                            )
            elif (l - elab) % 2 == 1:
                if not line.startswith(">"):
                    ilog.vlprint(
                        "Fasta file is either in Multi-line format (Single-line"
                        " format is required) or contains newlines or other "
                        "entries between sequence records.", vb-1
                        )
                    corrupt = True
                    if interact:
                        answer = input(
                            "Integration_locator can generate a converted copy"
                            " of the input file in the required format. Do you"
                            " wish to proceed (Y/N)?"
                            ).upper()
                    else:
                        ilog.vlprint(
                            "Integration_locator will generate a converted copy"
                            " of the input file in the required format.", vb
                        )
                        answer = 'Y'
                    if answer == 'N':
                        sys.exit()
                    elif answer == 'Y':
                        converted = True
                        try:
                            new_name = convert_multiline_fasta(_fasta)
                            corrupt = False
                        except:
                            ilog.vlprint(
                                "An unexpected error occurred while attempting "
                                "to generate converted input file. Please check"
                                "input fasta for irregularities", vb-1
                                )
                            raise
                        break
                    else:
                        ilog.vlprint('Invalid entry. Aborting program', vb-1)
                        sys.exit()
                else:
                    continue
            else:
                continue

    if fq_conv:
        ilog.vlprint(
            "Conducting extensive check of the converted file. should this"
            " check fail additional problems must have been present within your"
            " input fastq. In this case please carefully examine input files."
            , vb-1
            )
        corrupt, nc, nn, any_seq, any_def = check_input_file(
            new_name, "extensive", interact=False
            )

    if not converted:
        if elab > 0:
            ilog.vlprint("Attempting creation of fasta without empty lines", vb-1)
            corrupt = True
            converted = True
            try:
                new_name = convert_multiline_fasta(_fasta)
                corrupt = False
                ilog.vlprint(
                    "Successfully converted to file without empty lines. If"
                    " all fasta headers and sequence information were present"
                    " in the right order this will have produced a valid"
                    " file. However, if information was missing or disordered"
                    " in the source file the result will still be corrupt.", 
                    vb-1
                    )
            except:
                ilog.vlprint(
                    "ERROR: An unexpected error occurred while attempting "
                    "to generate converted input file. Please check"
                    "input fasta for irregularities", vb-1
                    )
                sys.exit()
        elif construct_fasta:
            new_name, converted = check_header_format(_fasta, l=1, h=1)
        else:
            pass
    
    if any_seq < 1:
        corrupt = True
        ilog.vlprint(f'File {_fasta} contains no sequence information', vb)
    if any_def < 1:
        corrupt = True
        ilog.vlprint(f'File {_fasta} file contains no fasta headers', vb)

    return corrupt, converted, new_name, any_seq, any_def


def test_input(_in_file, _check, constr_fa=False):
    """Test if input file is corrupt or requires reformatting"""
    check_input = check_input_file(_in_file, _check, construct_fasta=constr_fa)
    if check_input[0]:
        sys.exit()
    elif check_input[1]:
        _in_file = check_input[2]
    else:
        pass
    return _in_file


def chk_assembly_ver(species, ref_genome, args):
    """Automatic identification and matching of reference assembly"""
    ilog.vlprint(
        "Starting Automatic identification and matching of reference assembly",
        2
        )
    assembly_unknown = False
    method = "blastdbcmd"
    if not any([args.quick_mm2, args.intra]):
        try:
            check_for_refdb(ref_genome, args)
            refgen_info = get_genome_info_blast_tool(ref_genome)
        except:
            ilog.vlprint(
                f"INFO: get_genome_info_blast for {ref_genome} failed in "
                f"chk_assembly_ver. Attempting scan_ref_genome.", 1
            )
            refgen_info = scan_ref_genome_tool(ref_genome, {}, standalone=True)
            method = "scanref"
    else:
        ilog.vlprint(f"Using scan_ref_genome to get genome info.", 1)
        refgen_info = scan_ref_genome_tool(ref_genome, {}, standalone=True)
        method = "scanref"
    pkg_dir = os.path.split(__file__)[0]
    with open(os.path.join(pkg_dir, "il_resources", "genome_ident_db.json")) as inf:
        db_js = json.load(inf)
    spec_asms = db_js[species]
    categories = ["IDs", "len", "descr", "hash"]
    inv_len_dict = defaultdict(list)
    inv_hash_dict = defaultdict(list)
    inv_descr_dict = defaultdict(list)
    id_dict = defaultdict(list)
    for asm in spec_asms:
        id_dict[asm] = set(spec_asms[asm].keys())
        for tig in spec_asms[asm]:
            inv_len_dict[spec_asms[asm][tig]["len"][method]].append(asm)
            inv_descr_dict[spec_asms[asm][tig]["descr"][method]].append(asm)
            inv_hash_dict[spec_asms[asm][tig]["hash"][method]].append(asm)
    best_matches = {}
    perfect_matches = {}
    query_ids = set(refgen_info.keys())
    # compare ids
    best_match_ids = max(
        {
            asm: query_ids.intersection(id_dict[asm]) for asm in id_dict
            }.items(), key=lambda x: len(x[1])
        )
    if len(best_match_ids[1]) > 0:
        best_matches["IDs"] = ((best_match_ids[0], len(best_match_ids[1])),)
    id_match_ratio = len(best_match_ids[1])/len(query_ids)
    if id_match_ratio == 1:
        perfect_matches["IDs"] = best_match_ids
    
    def compare_category(catg, inv_cat_dct, refgen_info, spec_asms):
        query_cat = [refgen_info[tig][catg][method] for tig in refgen_info]
        asm_hit_ct = Counter({asm: 0 for asm in spec_asms})
        for l in query_cat:
            asm_hits = set(inv_cat_dct[l])
            asm_hit_ct.update(asm_hits)
        best_match_cat = max(asm_hit_ct.items(), key=lambda x: x[1])
        best_cat_ct = [
        (k, v) for k, v in asm_hit_ct.items() if v == best_match_cat[1]
        ]
        cat_match_ratio = best_match_cat[1]/len(query_ids)
        if len(best_cat_ct) > 1:
            above_thr = [
                match for match in best_cat_ct if match[1]/len(query_ids) > 0.1
                ]
            if any(above_thr):
                best_matches[catg] = tuple(best_cat_ct)
        else:
            # if best_match_cat[1] > 0:
            if cat_match_ratio > 0.1:
                best_matches[catg] = (best_match_cat,)
        
        
        if cat_match_ratio == 1:
            perfect_matches[catg] = best_match_cat[0]

    compare_category("len", inv_len_dict, refgen_info, spec_asms)
    compare_category("descr", inv_descr_dict, refgen_info, spec_asms)
    compare_category("hash", inv_hash_dict, refgen_info, spec_asms)

    single_hits = Counter()
    multiple_hits = {}
    for categ in best_matches:
        if len(best_matches[categ]) == 1:
            single_hits[best_matches[categ][0]] += 1
        elif len(best_matches[categ]) > 1:
            multiple_hits[categ] = best_matches[categ]
        else:
            ilog.vlprint("WARNING: Unexpected key type in chk_assembly_ver", 0)
    # for asm in single_hits:
    #     single_hits[cat] = tuple(
    #         [(asm,ct) for asm,ct in single_hits[cat] if ct > 0]
    #         )
    single_hits = {k:v for k,v in single_hits.items() if k[1] > 0}
    for cat in multiple_hits:
        multiple_hits[cat] = tuple(
            [(asm,ct) for asm,ct in multiple_hits[cat] if ct > 0]
            )
    multiple_hits = {k:v for k,v in multiple_hits.items() if len(v) > 0}

    if len(single_hits) > 0:
        best_single_hit = max(single_hits.items(), key=lambda x: x[1])[0][0]
        for cat in multiple_hits:
            for hit in multiple_hits[cat]:
                if hit[0] == best_single_hit:
                    best_matches[cat] = (hit,)
    elif len(single_hits) == 0 and len(multiple_hits) > 0:
        most_freq_hit = Counter()
        for cat in multiple_hits:
            asm_names = [hit[0] for hit in multiple_hits[cat]]
            most_freq_hit.update(asm_names)
        max_freq_hit = max(most_freq_hit.items(), key=lambda x: x[1])
        best_freq_hits = [
        (k, v) for k, v in most_freq_hit.items() if v == max_freq_hit[1]
        ]
        if len(best_freq_hits) == 1:
            hfq_hit = best_freq_hits[0][0]
        else:

            hfq_hit = max(most_freq_hit, key=lambda x: int(x.split(".")[-1]))
        for cat in multiple_hits:
            if hfq_hit in [hit[0] for hit in multiple_hits[cat]]:
                for hit in multiple_hits[cat]:
                    if hit[0] == hfq_hit:
                        best_matches[cat] = (hit,)
            else:
                best_matches[cat] = (multiple_hits[-1],)
    else:
        assembly_unknown = True
        ilog.vlprint(
            f"WARNING: The reference genome you are using could not be matched"
            f" to any assembly of {species} that intloc is aware of." 
            f" Please check if you indicated the right species and that you"
            f" chose the right reference file."
            f" A species specific ideogram can not be generated in this case."
            f" The relative positions of integrations in your reference will"
            f" still be reported. In order to enable intloc to create a high" 
            f" quality ideogram depiction it is recomended to retrieve an"
            f" assembly from NCBI (RefSeq, GenBank), Ensembl or a similary"
            f" renowned source.", 0
            )

    if assembly_unknown:
        return False, False
    else:
        # unpack outer tuple
        best_matches = {k: v[0] for k,v in best_matches.items()}
        for c in best_matches:
            ilog.vlprint(
                f"auto-ident assembly: best match {c} is {best_matches[c][0]} "
                f"with {best_matches[c][1]}/{len(query_ids)} matches", 2
                )
        for c in categories:
            if c not in best_matches:
                ilog.vlprint(
                    f"auto-ident assembly: no matches for category {c} "
                    f"were found", 3
                    )
        # select only asm name
        best_matches = {k: v[0] for k,v in best_matches.items()}

        best_overall_match = max(
            Counter(best_matches.values()).items(), key=lambda x: x[0]
            )

        perfect_match = True if len(perfect_matches) == 4 else False 

        
        ilog.vlprint(
            f"auto-ident assembly: The RefSeq assembly best matching the used "
            f"reference sequence is {best_overall_match[0]} with "
            f"{best_overall_match[1]}/{len(best_matches)} best matches ", 2
            )
        if perfect_match:
            ilog.vlprint(
                f"auto-ident assembly: The equivalent RefSeq assembly was "
                f"unequivocally identified by perfect match", 2
                )

        return best_overall_match, perfect_match


def check_header_format(_fasta, l=1, h=1):
    converted = False
    new_name = _fasta
    lc = 0
    hc = 0
    format_issue = False
    with open(_fasta, "r") as fin:
        for line in fin:
            lc += 1
            if lc > l:
                break
            elif hc > h:
                break
            elif line.startswith(">"):
                hc += 1
                new_header = fix_header(line)
                if new_header == line:
                    pass
                else:
                    format_issue = True
    
    if format_issue:
        ilog.vlprint(
            f"Fasta-headers in {_fasta} do not comply with required format."
            f" Generating reformatted filecopy for intloc."
            , 2
            )
        new_name = convert_multiline_fasta(_fasta)
        converted = True
    return new_name, converted


def fix_header(header):
    header = header.replace(", ", "_")
    header = header.replace(",", "_")
    seqname = header[1:].strip().split(" ")
    return f">{'_'.join(seqname)}\n"


def get_reads_by_header(search_list, all_reads_in, sel_read_out, specifier=""):
    """Retrieve reads from fasta file using list of headers"""
    ilog.vlprint(f'Getting {specifier} candidate reads using header info', 2)
    search_list = [r.strip() for r in search_list]
    chk_gt_no =  [r for r in search_list if not r.startswith(">")]
    chk_gt_yes =  [r for r in search_list if r.startswith(">")]
    if any(chk_gt_no) and not any(chk_gt_yes):
        search_list = [f">{r}" for r in search_list]
    elif any(chk_gt_no) and any(chk_gt_yes):
        ilog.vlprint("WARNING: Search_list contains mixed readname formats", 0)
    else:
        pass
    
    found_count = Counter({header: 0 for header in search_list})
    match = False
    int_reads = open(sel_read_out, "a")
    with open(all_reads_in, 'r') as fileo:
        for line in fileo:
            if match is False:
                if line.startswith('>'):
                    line_tmp = line.strip()
                    hls = line_tmp.split(" ")[0]
                    if hls in search_list:
                        found_count[hls] += 1
                        not_accepted = [",", ":", "|", ";", "\\", "/"]
                        lt_scan = [sym for sym in line_tmp if sym in not_accepted]
                        if len(lt_scan) > 0:
                            for sym in lt_scan:
                                line_tmp = line_tmp.replace(sym, " ")
                        line_tmp = line_tmp.strip().split(" ")
                        # print(line_tmp[0] + "\n", end='', file=int_reads)
                        int_reads.write(line_tmp[0] + "\n")
                        match = True
                        # testing if condensing search list speeds up
                        search_list.remove(hls)
                else:
                    pass
            elif match:
                # print(line, end='', file=int_reads)
                int_reads.write(line)
                match = False
    int_reads.close()

    return found_count


def convert_multiline_fasta(_fasta):
    """Converts multi-line fasta to single-line format"""
    ilog.vlprint('Converting multi-line fasta to single-line format', 3)

    head, tail = os.path.split(_fasta)
    spt = tail.split(".")
    new_name = "".join(["".join([i + "_" for i in spt[:-1]]), "sl.", spt[-1]])
    new_path = os.path.join(head, new_name)
    if os.path.exists(new_path):
        ilog.vlprint(
            "A converted single-line fasta copy appears to already exist.", 1)
    else:
        converted = open(new_path, "a")
        row = 0
        with open(_fasta, "r") as multiFile:
            for line in multiFile:
                if row == 0 and line.startswith(">"):
                    line = fix_header(line)
                    print(line, end="", file=converted)
                    row += 1
                elif row > 0 and line.startswith(">"):
                    line = fix_header(line)
                    print("\n" + line, end="", file=converted)
                elif len(line) == 0:
                    pass
                else:
                    print(line.strip(), end="", file=converted),
        converted.close()
        ilog.vlprint(
            "Input multi-line fasta read file was converted to single-line "
            "fasta.", 1
            )
    return os.path.join(head, new_name)


def sl_fasta_generator(fasta):
    with open(fasta, "r") as fa:
        for line in fa:
            if line.startswith(">"):
                yield (line.strip(), next(fa).strip())


def sl_fasta2dict(fasta):
    "load single-line fasta data into dict (from small files)"
    read_dict = dict()
    fa_gen = sl_fasta_generator(fasta)
    for id,seq in fa_gen:
        read_dict[id] = seq
    return read_dict


def fastq_to_fasta(fastq_in):
    """Converts fastq to fasta format"""
    ilog.vlprint('Converting fastq file to fasta format', 3)

    head, tail = os.path.split(fastq_in)
    spt = tail.split(".")
    new_name = "".join(["".join([i + "_" for i in spt[:-1]]), "conv.fasta"])

    def cont_gen():
        while True:
            yield "header"
            yield "seq"
            yield "plus"
            yield "qstring"

    content = cont_gen()
    bases = ["A", "T", "G", "C", "N"]
    with open(fastq_in, "r") as fastq:
        fasta = open(os.path.join(head, new_name), "w")
        for line in fastq:
            cur_line = next(content)
            if cur_line == "header":
                if line.startswith("@"):
                    line = line.replace("@", ">")
                    fasta.write(line)
            elif cur_line == "seq":
                if line[0].upper() in bases:
                    fasta.write(line)
                else:
                    ilog.vlprint("Unexpected character in sequence line", 3)
            elif cur_line == "plus":
                if line.startswith("+"):
                    pass
                else:
                    ilog.vlprint("Unexpected character in plus line",  3)
            else:
                pass
        fasta.close()

    return os.path.join(head, new_name)


if __name__ == "__main__":
    species = sys.argv[1]
    ref_genome = sys.argv[2]
    chk_assembly_ver(species, ref_genome)