# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import os
from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

def read_il_csv(csv_file, fdir=False):
    """read intloc csv files and return dict"""
    # file_dir path can be provided as string or list
    if fdir:
        if type(fdir)==list:
            fdir = os.path.join(*fdir)
        csv_file = os.path.join(fdir, csv_file)

    with open(csv_file, "r") as csv:
        data = csv.read().split("\n")
        categories, data = data[0].strip().split(","), data[1:]
        categories = [cat.strip() for cat in categories]
        data = [line for line in data if len(line) > 0]
        data = [
            {k:v for k,v in zip(categories, line.strip().split(","))} 
            for line in data
            ]
    return data


def load_data(out_dir, args):
    """Loading data for creation of html report"""
    ilog.vlprint("Loading data for creation of html report", 2)
    rep_data = {}

    run_file = os.path.join(out_dir, "run_info.csv")
    with open(run_file, "r") as runin:
        run_info = runin.read().split("\n")
        categories = ["category", "run-info"]
        run_info = [line.strip().split(",") for line in run_info]
        run_info = [
            {k:v for k,v in zip(categories, line)}
            for line in run_info if len(line)==2
            ]
        rep_data["run_info"] = run_info
    
    int_info = read_il_csv("Integration_Report.csv", fdir=out_dir)
    if not any([args.quick_mm2, args.polyclonal, args.intra]): # standard mode
        excl_info = ["read_names", "INT_name", "rs(m)"]
    elif args.quick_mm2:
        excl_info = ["read_names", "INT_name", "int-rs"]
    elif args.polyclonal:
        excl_info = ["read_names", "INT_name", "int-rs"]
    elif args.intra:
        excl_info = [
            "read_names",
            "INT_name",
            "int-rs",
            "+/-[bp]",
            "non-int-r",
            "ts_cov",
            "IMR",
            "SSBR"
            ]
    else:
        excl_info = ["read_names", "INT_name"]
    
    for dct in int_info:
        for info in excl_info:
            dct.pop(info)
    rep_data["int_info"] = int_info

    args.files = [e.name for e in os.scandir()]
    if args.species:
        feat_info = read_il_csv("int_features.csv", fdir=out_dir)
        rep_data["feat_info"] = feat_info
    # if not any(
    #     [
    #         args.no_multi_ints,
    #         args.skip_eval,
    #         args.quick_mm2,
    #         args.intra,
    #         args.polyclonal
    #         ]
    #         ):
    if "proximal_trans_sites_report.csv" in args.files:
        proxtrans_info = read_il_csv("proximal_trans_sites_report.csv", fdir=out_dir)
        rep_data["proxtrans_info"] = proxtrans_info
        multisite_info = read_il_csv("multi_site_report.csv", fdir=out_dir)
        rep_data["multisite_info"] = multisite_info
    if args.intra:
        preex_info = read_il_csv("preex_Integration_Report.csv", fdir=out_dir)
        for dct in preex_info:
            # dct.pop("read_names")
            for info in excl_info:
                dct.pop(info)
        rep_data["preex_int_info"] = preex_info
        if args.species:
            preex_feat_info = read_il_csv("preex_int_features.csv", fdir=out_dir)
            rep_data["preex_feat_info"] = preex_feat_info
    if args.cluster_reads:
        clust_info = read_il_csv("ava_read_cluster.csv", fdir=out_dir)
        for dct in clust_info:
            dct["read names"] = ", ".join(dct["read names"].split(";"))
        rep_data["clust_info"] = clust_info
    if "seq_gain_loss_at_intsites.csv" in args.files:
        lossgain = read_il_csv("seq_gain_loss_at_intsites.csv", fdir=out_dir)
        rep_data["lossgain_info"] = lossgain

    return rep_data


def load_figures(out_dir, args):
    """Loading figures for creation of html report"""
    ilog.vlprint("Loading figures for creation of html report", 2)
    figures = {}

    if args.species:
        file = os.path.join(out_dir, "ideogram.svg")
        with open(file, "r") as idin:
            ideo = idin.read().split("\n")
    else:
        file = os.path.join(out_dir, "relative_integration_locations.svg")
        with open(file, "r") as idin:
            ideo = idin.read().split("\n")
    figures["ideo"] = ideo

    if args.intra:
        if args.species:
            file = os.path.join(out_dir, "preex_ideogram.svg")
            with open(file, "r") as idin:
                px_ideo = idin.read().split("\n")
        else:
            file = os.path.join(out_dir, "preex_relative_integration_locations.svg")
            with open(file, "r") as idin:
                px_ideo = idin.read().split("\n")
        figures["preex_ideo"] = px_ideo

    if args.cluster_reads:
        file = os.path.join(out_dir, "read_cluster_il-ava-clust.svg")
        with open(file, "r") as clustin:
            clust = clustin.read().split("\n")
        figures["clust"] = clust
    
    return figures


def insert_content(line, inserts, inc_ct=0, increment=0, pos=0, consec=False):
            insertions = line.count("@")/2
            ls = line.split("@")
            replace = {}
            if insertions > 0:
                for seg in ls:
                    if seg.startswith("~") or seg.startswith("#"):
                        quotes = True if seg.startswith("~") else False
                        idx = ls.index(seg)
                        replace[idx] = inserts[seg[1:]]
                for idx in replace:
                    ls.pop(idx)
                    if quotes:
                        if not consec:
                            try:
                                pos = float(replace[idx]) + increment
                            except ValueError:
                                pos = replace[idx]
                        else:
                            pos = pos + increment
                        ls.insert(idx, f"\"{pos}\"")
                    else:
                        ls.insert(idx, f"{replace[idx]}")
            new_line = "".join(ls)
            return new_line, inc_ct, pos


def write_block(report, block):
    for line in block:
        report.write(line)


def insert_read_info(block, run_info):
    count_lines = {k:v for k,v in enumerate(block)}
    run_info = {dct["category"]:dct["run-info"] for dct in run_info}
    for line in count_lines:
        nl_ict_pos = insert_content(count_lines[line], run_info)
        count_lines[line] = nl_ict_pos[0]
    mod_block = count_lines.values()
    return mod_block


def insert_fig_svg(block, svg):
    for line in block:
        if line.strip().startswith("@"):
            idx = block.index(line)
            ls = line.split("@")
            if ls[2] == "all":
                block.pop(idx)
                # add correct indentation level
                svg = [ls[0] + line for line in svg]
                block.insert(idx, "\n".join(svg))
    return block


def insert_table_data(rep, table_name, data, html_blox, select=False, path=""):
    """
    Function for insertion of data tables in html report. Used for gener
    ation of regular tables or multi-figure display boxes with drop-down
        selection using tables as frame. If select=False the data is insert
    ed from a dict generated in load_data. If select=True a subfolder (o
    r sub-sub etc.) can be specified via the path kwarg. All figure-file
    s in the subfolder will be included in the drop-down.
    """
    tbl_head = f"{table_name}_head"
    tbl_data = f"{table_name}_data"
    tbl_foot = f"{table_name}_foot"
    write_block(rep, html_blox[tbl_head])
    data_sl = html_blox[tbl_data]

    if not select:
        ind = html_blox[tbl_data][1].split("@")[0]
        for row_dct in data:
            rep.write(data_sl[0])
            for k in row_dct:
                rep.write(f"{ind}<td>{row_dct[k]}</td>\n")
            rep.write(data_sl[-1])
        write_block(rep, html_blox[tbl_foot])
    else:
        ind = html_blox[tbl_data][0].split("@")[0]
        fpath = ""
        for scanned_file in data:
            fname = scanned_file.name
            fpath = (os.path.join(".", path, fname))
            rep.write(ind + f"<option value=\"{fpath}\">{fname}</option>\n")
        footer = html_blox[tbl_foot]
        lc = -1
        for line in footer:
            lc += 1
            if "@" in line:
                ins_slot = lc
        sl = footer[ins_slot].split("@")
        sl.insert(1, fpath)
        footer.pop(ins_slot)
        footer.insert(ins_slot,"".join(sl))
        write_block(rep, footer)


def write_html_report(args):
    """Write html report"""
    ilog.vlprint("Writing html report", 1)
    outfile = os.path.join(args.outdir, "Integration_summary_report.html")
    pkg_dir = os.path.split(__file__)[0]
    template = os.path.join(pkg_dir, "il_resources", "html_rep_template.html")
    html_blox = {}
    files = [e.name for e in os.scandir()]

    with open(outfile, "w") as rep, open(template, "r") as templ:
        first_block = True
        lines = []
        for line in templ:
            if line.startswith("$"):
                if first_block:
                    first_block = False
                    key = line.strip().split("$")[1]
                else:
                    html_blox[key] = lines
                    key = line.strip().split("$")[1]
                    lines = []
            else:
                if first_block:
                    pass
                else:
                    lines.append(line)
    
        data = load_data(args.outdir, args)
        figures = load_figures(args.outdir, args)

                    
        write_block(rep, html_blox["preamble"])
        # mode_info = insert_read_info(html_blox["mode"], data["run_info"])
        # write_block(rep, mode_info)
        if args.intra:
            read_info = insert_read_info(html_blox["preex_read_info"], data["run_info"])
        else:
            read_info = insert_read_info(html_blox["read_info"], data["run_info"])
        write_block(rep, read_info)
        if args.species:
            if args.intra:
                ideo = insert_fig_svg(html_blox["novel_spec_ideo"], figures["ideo"])
            else:
                ideo = insert_fig_svg(html_blox["spec_ideo"], figures["ideo"])
        else:
            ideo = insert_fig_svg(html_blox["generic_ideo"], figures["ideo"])
        write_block(rep, ideo)
        if args.intra:
            insert_table_data(rep, "new_intras_table", data["int_info"], html_blox)
        else:
            insert_table_data(rep, "int_table", data["int_info"], html_blox)
        if "seq_gain_loss_at_intsites.csv" in files: 
            # add table with seq loss/gain at intsite
            insert_table_data(rep, "lossgain_table", data["lossgain_info"], html_blox)
        if args.species:
            singles = os.scandir("single_chrom")
            lf = os.path.join("figures", "single_chrom")
            insert_table_data(rep, "single_ideo", singles, html_blox, select=True, path=lf)
            insert_table_data(rep, "feat_table", data["feat_info"], html_blox)
        try:
            blast_maps = os.scandir("blast_map_figs")
            lf = os.path.join("figures", "blast_map_figs")
            insert_table_data(rep, "cand_read_maps_blast", blast_maps, html_blox, select=True, path=lf)
        except Exception as e:
            ilog.vlprint(e, 2)
            ilog.vlprint(
                "WARNING: Integration of blast candidate-read mapping figures"
                " into html report failed.", 2
                )
        if not args.skip_eval:
            # add optional minimap2 mapped reads figure
            try:
                mm2_maps = os.scandir("minimap2_map_figs")
                lf = os.path.join("figures", "minimap2_map_figs")
                insert_table_data(rep, "cand_read_maps_mm2", mm2_maps, html_blox, select=True, path=lf)
            except Exception as e:
                ilog.vlprint(e, 2)
                ilog.vlprint(
                    "WARNING: Integration of mm2 candidate-read mapping figures"
                    " into html report failed. Possibly due to failure of"
                    " il_evaluate.", 2
                    )
        # if all([not args.no_multi_ints, not args.skip_eval]):
        if "proxtrans_info" in data:
            if len(data["proxtrans_info"]) > 0:
                insert_table_data(rep, "proxtrans_table", data["proxtrans_info"], html_blox)
            if len(data["multisite_info"]) > 0:
                insert_table_data(rep, "multisite_table", data["multisite_info"], html_blox)
        if args.intra:
            if args.species:
                px_ideo = insert_fig_svg(html_blox["preex_spec_ideo"], figures["preex_ideo"])
                write_block(rep, px_ideo)
                insert_table_data(rep, "preex_int_table", data["preex_int_info"], html_blox)
                insert_table_data(rep, "preex_feat_table", data["preex_feat_info"], html_blox)
            else:
                px_ideo = insert_fig_svg(html_blox["preex_generic_ideo"], figures["preex_ideo"])
                write_block(rep, px_ideo)
                insert_table_data(rep, "preex_int_table", data["preex_int_info"], html_blox)
        if args.cluster_reads:
            clust_map = insert_fig_svg(html_blox["clust_fig"], figures["clust"])
            write_block(rep, clust_map)
            insert_table_data(rep, "clust_table", data["clust_info"], html_blox)
        insert_table_data(rep, "info_table", data["run_info"], html_blox)
        write_block(rep, html_blox["html_foot"])