import os
import sys
import subprocess
from il_logging import Ilogger

ilog = Ilogger()
ilog.module = __name__

def caplog_sp_error(error, message, exit=True):
    ilog.vlprint(f"Error: {message}", -1)
    print("Check intlog.log for error details")
    ilog.vlprint(error.stdout.decode('utf-8'), 2)
    ilog.vlprint(error.stderr.decode('utf-8'), 2)
    if exit:
        ilog.vlprint("Exiting integration locator", 1)
        sys.exit()


def minimap2_sp(
        target,
        query,
        mm_out,
        il_args,
        preset="-x",
        map_mode="map",
        secondary="yes",
        N="5",
        optargs=[],
        check=True,
        capture_output=True,
        ):
    command = [
        il_args.minimap2,
        f"--secondary={secondary}",
        "-N", N,
        preset, f"{map_mode}-{il_args.seq_tec}",
        target,
        query
    ]
    command.extend(optargs)

    if il_args.dep_ver["minimap2"]!="outdated":
        command.extend(["-o", mm_out])
        subprocess.run(command, check=check, capture_output=capture_output)
    else:
        mm_res = subprocess.run(
            command,
            check=check,
            capture_output=capture_output,
            )
        with open(mm_out, "w") as out:
            out.write(mm_res.stdout.decode())


def intloc_self_sp(cmd, outdir=False):
    """Spawn a subprocess of intloc itself"""
    intloc_path = os.path.split(__file__)[0]
    if not outdir:
        outdir = os.path.join(intloc_path, "il_test", "Test_results")
    compl_proc_obj = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return compl_proc_obj


# ####################il_intra######################################
# subprocess.run(
#             [
#             args.minimap2,
#             "-x", f"map-{args.seq_tec}",
#             "--secondary=no", # available in minimap2 2.3 and higher
#             args.construct,
#             genome,
#             "-o", intra_sites  # available in minimap2 2.16 and higher
#             ],
#             check=True, capture_output=True
#         )
# ############################################################################
# ############################################################################


# ####################il_polyclonal######################################
# subprocess.run(
#                 [
#                 args.minimap2,
#                 "-N", "0",
#                 "-x", f"map-{args.seq_tec}",
#                 genome,
#                 int_cand_reads,
#                 "-o", outfile
#                 ],
#             check=True, capture_output=True
#             )
# ############################################################################
# ############################################################################

# ####################integration_locator(main)###################################
# subprocess.run(
#             [
#             args.blastn,
#             "-db", genome,
#             "-query", construct,
#             "-out", tmp,
#             "-outfmt", "10 qaccver qlen saccver slen length qstart qend sstart"
#             " send pident mismatch gapopen evalue bitscore qcovhsp qcovus",
#             "-num_threads", args.cores,
#             "-culling_limit", "1",
#             "-evalue", "0.0001",
#             "-perc_identity", "95"
#             ],
#             check=True, capture_output=True
#             )


# subprocess.run(
#                 [
#                 args.makeblastdb,
#                 "-in", reads,
#                 "-parse_seqids",
#                 "-dbtype", "nucl"
#                 ],
#                 check=True, capture_output=True
#                 )


# subprocess.run([
#                 args.blastn,
#                 "-db", reads,
#                 "-query", construct,
#                 "-out", tmp,
#                 "-outfmt", 
#                 "10 qaccver qlen saccver slen length qstart qend sstart send "
#                 "pident mismatch gapopen evalue bitscore qcovhsp qcovus sstrand",
#                 "-num_threads", args.cores,
#                 "-evalue", "0.000000000000000001",
#                 "-perc_identity", "80",
#                 "-dust", "no",
#                 "-max_target_seqs", "10000000"],
#                 check=True, capture_output=True)

# subprocess.run([
#                 args.blastn,
#                 "-db", genome,
#                 "-query", int_cand_reads,
#                 "-out", tmp,
#                 "-outfmt", "10 qaccver qlen saccver slen length "
#                 "qstart qend sstart send pident mismatch gapopen "
#                 "evalue bitscore qcovhsp qcovus score gaps nident "
#                 "positive sstrand", 
#                 "-num_threads", args.cores,
#                 "-culling_limit", "1",
#                 "-evalue", "0.0001",
#                 "-perc_identity", "80",
#                 "-gapopen", "0",
#                 "-gapextend", "0"],
#                 check=True, capture_output=True
#                 )
# ############################################################################
# ############################################################################    