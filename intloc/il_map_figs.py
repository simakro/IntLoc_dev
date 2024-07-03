import os

from il_html import insert_content
from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

def gmf_write_block(figure, block):
        for line in block:
            figure.write(line)


def gmf_load_block(block):
    class_lines = {}
    count_lines = {}
    lc = 0
    for line in block:
        lc += 1
        count_lines[lc] = line
        if "class=" in line:
            ls_class = line.split("class=\"")
            key = ls_class[1].split("\"")[0]
            class_lines[key] = line
    return class_lines, count_lines


# def gmf_write_preamble(int_site, preamble, aligner="blast"):
#     class_lines, count_lines = gmf_load_block(preamble)
#     span_reg = int_site.spanned_region[aligner]
#     al_ct = len(int_site.supp_als) if aligner=="blast" else len(int_site.minimap_supp_sr)
#     al_ct = al_ct if al_ct<200 else 200
#     inserts = {}
#     inserts["ref_len"] = (span_reg[1] - span_reg[0] + 1) / 100
#     inserts["svg_height"] = al_ct * 15 + 60
#     inserts["svg_width"] = inserts["ref_len"] + 40
#     for line in count_lines:
#         nl_ict_pos = insert_content(count_lines[line], inserts)
#         count_lines[line] = nl_ict_pos[0]
#     count_lines = list(count_lines.values())
#     return count_lines


def gmf_write_preamble(int_site, preamble, aligner="blast", fig="map_fig"):
    class_lines, count_lines = gmf_load_block(preamble)
    inserts = {}
    if fig=="map_fig":
        span_reg = int_site.spanned_region[aligner]
        al_ct = len(int_site.supp_als) if aligner=="blast" else len(int_site.minimap_supp_sr)
        al_ct = al_ct if al_ct<200 else 200
        inserts["ref_len"] = (span_reg[1] - span_reg[0] + 1) / 100
        inserts["svg_height"] = al_ct * 15 + 60
        inserts["svg_width"] = inserts["ref_len"] + 40
    else:
        read_ct = len(int_site)
        max_len = max([cro.slen for cro in int_site.values()])/100
        inserts["svg_height"] = read_ct * 25 + 60
        inserts["svg_width"] = max_len + 40
    for line in count_lines:
        nl_ict_pos = insert_content(count_lines[line], inserts)
        count_lines[line] = nl_ict_pos[0]
    count_lines = list(count_lines.values())
    return count_lines


def gmf_write_reference(int_site, ref_block, aligner="blast"):
    class_lines, count_lines = gmf_load_block(ref_block)
    span_reg = int_site.spanned_region[aligner]
    inserts = {}
    start_span = span_reg[0]
    end_span = span_reg[1]
    inserts["max_start_label"] = start_span
    inserts["pos_label"] = 10
    inserts["ref_len"] = str((span_reg[1] - span_reg[0] + 1) / 100)
    inserts["int_mark"] = (int_site.approx_loc-start_span)/100 + 20
    inserts["int_label"] = (int_site.approx_loc-start_span)/100 - 5
    inserts["major_x"] = 20
    inserts["minor_x"] = 40

    for line in count_lines:
        nl_ict_pos = insert_content(count_lines[line], inserts)
        count_lines[line] = nl_ict_pos[0]

    base_line = list(count_lines.values())
    last_pos = 40
    ict = 3
    while last_pos < (end_span-start_span)/100-20:
        if ict % 5 == 1:
            ict += 1
            nl,ict, pos = insert_content(
                class_lines["major_tick"], inserts, 
                inc_ct=ict, increment=20, consec=True, pos=last_pos
                )
            inserts["pos_label"] += 100
            inserts["max_start_label"] += 10000

            mtick_label = insert_content(
                class_lines["tick_txt"], inserts, 
                inc_ct=ict, increment=20, consec=True, pos=last_pos
                )
            nl = nl + mtick_label[0]
        else:
            ict += 1
            nl,ict, pos = insert_content(
                class_lines["minor_tick"], inserts, inc_ct=ict, increment=20, consec=True, pos=last_pos
                )
        last_pos = pos
        base_line.insert(-2, nl)
    return base_line


def gmf_write_reads(int_site, ref_block, aligner="minimap2"):
    class_lines, count_lines = gmf_load_block(ref_block)
    span_reg = int_site.spanned_region[aligner]
    inserts = {}
    start_span = span_reg[0]
    read_ct = 0
    read_z = 30

    # if aligner =="minimap2":
    #     container = int_site.minimap_supp_sr
    # elif aligner =="blast":
    #     container = int_site.supp_als
    # else:
    #     # this is reserved for multi-int-read plotting
    #     container = int_site
    container = int_site.minimap_supp_sr if aligner=="minimap2" else int_site.supp_als
    read_paragraphs = []
    for read in container:
        read_ct += 1
        split_ct = 0
        coord_source = container[read] if aligner=="minimap2" else read.split_coord
        for aln in coord_source:
            split_ct += 1
            paragraph = []
            readname = aln["qname"] if aligner=="minimap2" else read.READ
            tstart = aln["tstart"]+start_span if aligner=="minimap2" else min(aln) 
            tend = aln["tend"]+start_span if aligner=="minimap2" else max(aln) 
            inserts["al_len"] = aln["map_len"]/100 if aligner=="minimap2" else (abs(tend-tstart)+1)/100
            inserts["ttip_txt"] = f"{readname}\\ntstart:{tstart} tend:{tend}"
            inserts["read_x"] = (aln["tstart"])/100+20 if aligner=="minimap2" else (tstart-start_span)/100+20
            inserts["read_z"] = read_z + read_ct * 15
            inserts["color"] = "#008000" if split_ct == 1 else "#00a000"
            minus = "-" if aligner=="minimap2" else "minus"
            strand = aln["strand"] if aligner=="minimap2" else read.sstrand[0]
            direct = 0 if strand==minus else inserts["al_len"]
            inserts["direction"] = direct
            arrow = "5 0, 0 5, 5 10" if direct==0 else "-5 0, 0 5, -5 10"
            inserts["arrowhead"] = arrow
            for line in count_lines:
                nl_ict_pos = insert_content(count_lines[line], inserts)
                paragraph.append(nl_ict_pos[0])
            read_paragraphs.extend(paragraph)
    return read_paragraphs


def gen_map_figs(primary_int_site_objs, aligner="blast"):
    """Draw figures depicting mappings of reads to reference with aligner."""
    ilog.vlprint(f'Drawing figures of read alignments to ref with {aligner}', 3)

    pkg_dir = os.path.split(__file__)[0]
    template = os.path.join(pkg_dir, "il_resources", "templ_read_map_fig.html")
    svg_blox = {}

    with open(template, "r") as templ:
        first_block = True
        lines = []
        for line in templ:
            if line.startswith("$"):
                if first_block:
                    first_block = False
                    key = line.strip().split("$")[1]
                else:
                    svg_blox[key] = lines
                    key = line.strip().split("$")[1]
                    lines = []
            else:
                lines.append(line)
        svg_blox[key] = lines

    try:
        os.mkdir(f"{aligner}_map_figs")
    except FileExistsError:
        ilog.vlprint(f"Directory {aligner}_map_figs already exists", 1)
    for piso in primary_int_site_objs:
        outf = os.path.join(f"{aligner}_map_figs", piso.name + "_cand_read_mappings.html")
        try:
            with open(outf, "w") as out:
                preamble = gmf_write_preamble(piso, svg_blox["preamble"], aligner=aligner)
                gmf_write_block(out, preamble)
                reference = gmf_write_reference(piso, svg_blox["reference"], aligner=aligner)
                gmf_write_block(out, reference)
                int_reads = gmf_write_reads(piso, svg_blox["int_read"], aligner=aligner)
                gmf_write_block(out, int_reads)
                # write_block(out, svg_blox["non_int_read"])
                gmf_write_block(out, svg_blox["end_of_file"])
        except:
            ilog.vlprint(f"{aligner} read map figure creation for {piso.name} failed", 1)



def gmirf_write_multi_reads(mult_cro_dct, read_block, int_block, step="read"):
    count_lines = gmf_load_block(read_block)[1]
    cl_int_block = gmf_load_block(int_block)[1]
    inserts = {}
    read_ct = 0
    read_z = 30

    read_paragraphs = []
    for read in mult_cro_dct:
        read = mult_cro_dct[read]
        read_ct += 1
        # split_ct = 0
        # coord_source = container[read] if aligner=="minimap2" else read.split_coord
        # for aln in coord_source:
        # split_ct += 1
        paragraph = []
        readname = read.RNAME # if step=="read" else "Integration sequence"
        length = read.slen
        inserts["read_name"] = readname
        inserts["al_strand"] = "no al." if not read.rtr else read.rtr.sstrand[0]
        inserts["read_len"] = length/100
        inserts["ttip_txt"] = f"{readname}\\nlength:{length}"
        inserts["read_x"] = 20
        inserts["read_z"] = read_z + read_ct * 25
        inserts["color"] = "#00a000"
        direct = 5
        inserts["direction"] = direct
        arrow = ""
        inserts["arrowhead"] = arrow
        for line in count_lines:
            nl_ict_pos = insert_content(count_lines[line], inserts)
            paragraph.append(nl_ict_pos[0])
        read_paragraphs.extend(paragraph)
        for int_seq in read.orig_alignment:
            paragraph = []
            readname = "Integration sequence"
            rstart, rend = int_seq
            length = rend-rstart+1
            inserts["read_len"] = length/100
            inserts["ttip_txt"] = f"{readname}\\npos. in read:{rstart}-{rend}"
            inserts["read_x"] = 20 + rstart/100
            inserts["read_z"] = read_z + read_ct * 25
            inserts["color"] =  "#ff0000"
            # minus = "-" if aligner=="minimap2" else "minus"
            # strand = aln["strand"] if aligner=="minimap2" else read.sstrand[0]
            # direct = 0 if strand==minus else inserts["al_len"]
            try:
                strand = read.strand[tuple(int_seq)]
            except KeyError as ke:
                strand = "n.a."
                ilog.vlprint(
                    f"{ke} Alignment direction for int-seq {int_seq} of read"
                    f" {read.RNAME} could not be determined" ,2
                    )
            direct = 0 if strand=="minus" else (rend-rstart)/100
            inserts["direction"] = direct
            arrow = "5 0, 0 5, 5 10" if direct==0 else "-5 0, 0 5, -5 10"
            inserts["arrowhead"] = arrow
            for line in cl_int_block:
                nl_ict_pos = insert_content(cl_int_block[line], inserts)
                paragraph.append(nl_ict_pos[0])
            read_paragraphs.extend(paragraph)

    return read_paragraphs


def gen_multi_int_read_figs(mult_cro_dct, aligner="blast"):
    """Draw figures depicting position of integrating sequences in reads."""
    ilog.vlprint(
        f'Drawing figs depicting position of integrating sequences in reads', 3
        )

    pkg_dir = os.path.split(__file__)[0]
    template = os.path.join(pkg_dir, "il_resources", "templ_mutli_int_read.html")
    svg_blox = {}

    with open(template, "r") as templ:
        first_block = True
        lines = []
        for line in templ:
            if line.startswith("$"):
                if first_block:
                    first_block = False
                    key = line.strip().split("$")[1]
                else:
                    svg_blox[key] = lines
                    key = line.strip().split("$")[1]
                    lines = []
            else:
                lines.append(line)
        svg_blox[key] = lines

    outf = "multi_intread_fig.html"
    # try:
    with open(outf, "w") as out:
        preamble = gmf_write_preamble(mult_cro_dct, svg_blox["preamble"], aligner=aligner, fig="mult_int_reads")
        gmf_write_block(out, preamble)

        # reference = gmf_write_reference(piso, svg_blox["reference"], aligner=aligner)
        # gmf_write_block(out, reference)
        # int_reads = gmf_write_reads(piso, svg_blox["read"], aligner=aligner)
        int_reads = gmirf_write_multi_reads(mult_cro_dct, svg_blox["read"], svg_blox["int_seq"])
        gmf_write_block(out, int_reads)
        
        # int_seqs = gmf_write_reads(piso, svg_blox["int_seq"], aligner=aligner)
        # gmf_write_block(out, int_seqs)


        # write_block(out, svg_blox["non_int_read"])
        gmf_write_block(out, svg_blox["end_of_file"])
    # except:
    #     ilog.vlprint(f"Multi-int-read figure creation failed", 1)

