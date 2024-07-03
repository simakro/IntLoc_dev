# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import statistics as stat
import os
import json
import math
from collections import defaultdict

from il_logging import Ilogger
from il_chk_conv import report_time


ilog = Ilogger()
ilog.module = __name__

class Chromosome:
    def __init__(
        self,
        name,
        _block_pos,
        p_arm_svg,
        q_arm_svg,
        len_bp_total,
        template,
        chr_type,
        int_height,
        txt_displ,
        single_fig,
        args
        ):
        self.name = name
        self.block_pos = _block_pos
        self.p_arm_svg = p_arm_svg
        self.q_arm_svg = q_arm_svg
        if chr_type == "two_arms":
            self.total_len_svg = self.p_arm_svg + self.q_arm_svg
        else:
            self.total_len_svg = self.p_arm_svg
        self.len_bp_total = len_bp_total
        self.template = template
        self.int_height = int_height
        self.bp_svg_ratio = self.len_bp_total / self.total_len_svg
        self.block = self.template[self.block_pos[0]: self.block_pos[1]+1]
        self.integrations = []
        self.int_info = {}
        self.max_chrom_len = int()
        # attributes for localization of hover-over display text#
        # mcf=matrix_correction_factor, cs=chromosome_spacing, ls=line_spacing
        self.txt_x = txt_displ["x"] - self.name * txt_displ["cs"]
        self.txt_y = txt_displ["y"]
        self.font_size = txt_displ["fs"]
        self.lin_space = txt_displ["ls"]
        self.mcf = txt_displ["mcf"]
        if txt_displ["aln"] == "top":
            self.bottom_aln = False
            self.transl = 0
        else:
            self.bottom_aln = True
            self.transl = self.total_len_svg
        self.single_fig = single_fig
        self.args = args


    def insert_int(self, line_spacer):
        int_color = self.args.color_ints
        ints_conv = dict(
            sorted(
                [(intloc, intloc/self.bp_svg_ratio) for intloc in self.integrations]
                )
        )
        ignore = ["arm", "</g", "<g", "<text", "<use"]
        lc = 0
         
        for intg in ints_conv:
            info = self.int_info[intg]
            if self.args.intra:
                if int(info["mms"]) < 1:
                    int_color = "blue"
                else:
                    int_color = "lime"
            open_group = "              <g>\n"
            int_line = f'                <rect fill="{int_color}" class="Integration"'\
                f' height="{self.int_height}" stroke="none" width="12" x="0"'\
                f' y="{ints_conv[intg]}"/>\n'
            # int_line = f'                <line stroke="{int_color}" class="Integration"'\
            #     f' stroke-width="{self.int_height}" x1="0" x2="12"'\
            #     f' y1="{ints_conv[intg]}" y2="{ints_conv[intg]}"/>\n'
            int_txt_1 = f'                <text class="int_txt" fill="black"'\
                f' text-decoration="underline" font-size="{self.font_size}px"'\
                f' transform="translate(0,{self.transl}),matrix(1 0 0 {self.mcf} 0 0)"'\
                f' x="{self.txt_x}" y="{self.txt_y}"'\
                f' font-weight="bold" >Int.{info["int_id"]}</text>\n'
            int_txt_2 = f'                <text class="int_txt" fill="black"'\
                f' font-size="{self.font_size}px"'\
                f' transform="translate(0,{self.transl}),matrix(1 0 0 {self.mcf} 0 0)"'\
                f' x="{self.txt_x}" y="{self.txt_y + self.lin_space}"'\
                f' >Name: {info["int_name"] }</text>\n'
            int_txt_3 = f'                <text class="int_txt" fill="black"'\
                f' font-size="{self.font_size}px"'\
                f' transform="translate(0,{self.transl}),matrix(1 0 0 {self.mcf} 0 0)"'\
                f' x="{self.txt_x}" y="{self.txt_y + 2*self.lin_space}"'\
                f' >Chrom: {info["chr_name"]}</text>\n'
            int_txt_4 = f'                <text class="int_txt" fill="black"'\
                f' font-size="{self.font_size}px"'\
                f' transform="translate(0,{self.transl}),matrix(1 0 0 {self.mcf} 0 0)"'\
                f' x="{self.txt_x}" y="{self.txt_y + 3*self.lin_space}"'\
                f' >Contig: {info["contig_name"]}</text>\n'
            int_txt_5 = f'                <text class="int_txt" fill="black"'\
                f' font-size="{self.font_size}px"'\
                f' transform="translate(0,{self.transl}),matrix(1 0 0 {self.mcf} 0 0)"'\
                f' x="{self.txt_x}" y="{self.txt_y + 4*self.lin_space}"'\
                f' >Position [bp]: {info["int_pos_bp"]}</text>\n'
            int_txt_6 = f'                <text class="int_txt" fill="black"'\
                f' font-size="{self.font_size}px"'\
                f' transform="translate(0,{self.transl}),matrix(1 0 0 {self.mcf} 0 0)"'\
                f' x="{self.txt_x}" y="{self.txt_y + 5*self.lin_space}"'\
                f' >Features: {info["features"]}</text>\n'
            close_group = "              </g>\n"
            y_pos = float(int_line.split("\"")[-2]) + self.int_height
            idx = -1
            clip_path_section = False
            int_section = False
            group_open = False
            for line in self.block:
                idx += 1
                ls = line.split("\"")
                if clip_path_section:
                    if "</g" in ls[0]:
                        clip_path_section = False
                        int_section = True
                elif int_section:
                    if "<g>" in ls[0]:
                        group_open = True
                    elif "</g" in ls[0]:
                        if group_open:
                            group_open = False
                        else:
                            self.block.insert(idx, open_group)
                            self.block.insert(idx + 1, int_line)
                            self.block.insert(idx + 2, int_txt_1)
                            self.block.insert(idx + 3, int_txt_2)
                            self.block.insert(idx + 4, int_txt_3)
                            self.block.insert(idx + 5, int_txt_4)
                            self.block.insert(idx + 6, int_txt_5)
                            self.block.insert(idx + 7, int_txt_6)
                            self.block.insert(idx + 8, close_group)
                            break
                    else:
                        pass
                else:
                    if "clip-path=" in ls[0]:
                        clip_path_section = True
                        # if float(ls[-2]) > y_pos:
                        #     self.block.insert(idx, int_line)
                        #     break
    

class SpeciesResources:
    def __init__(self, species, spec_dict, res_dir):
        self.spec_name = species
        self.svg_src = spec_dict["svg_src"]
        self.len_bp_real = spec_dict["len_bp_real"]
        self.gene_dens = spec_dict["gene_dens_params"]
        self.feat_src = os.path.join(res_dir, spec_dict["feat_src"])
        self.chr_tr_dict = spec_dict["chr_tr_dict"]
        self.corr_val = spec_dict["corr_val"]
        self.txt_displ = spec_dict["txt_displ"]
        self.single_fig = spec_dict["single_fig"]
        self.chrom_scale = None
        self.int_height = None


    def calc_int_height(self):
        self.int_height = self.chrom_scale / 100


def generate_spec_resource(species, res_dir):
    ilog.vlprint("Generating species resource", 3)
    with open(os.path.join(res_dir, "il_spec_resources.json"), "r") as infile:
        json_in = json.load(infile)
    spec_dict = json_in[species]

    spec_resource = SpeciesResources(
                                species,
                                spec_dict,
                                res_dir
                                )

    return spec_resource


def parse_svg_src(_svg_src, spec_res):
    ilog.vlprint("Parsing svg source", 3)
    with open(_svg_src, "r") as sin:
        templ = sin.readlines()

    block_pos = {}
    arm_len = {}
    lc = 0
    cc = 0
    for line in templ:
        lc += 1
        if "use" in line:
            if len(block_pos) == 0:
                cc += 1
                block_pos[cc] = [lc-2]
            else:
                cc += 1
                block_pos[cc] = [lc-2]
                block_pos[cc-1].append(lc-3)
        elif line == '</svg>\n':
            block_pos[cc].append(lc)
        elif "arm" in line:
            if "linearGradient" in line:
                pass
            else:
                if cc in arm_len:
                    ls = line.split("\"")
                    arm_len[cc].append(float(ls[3]))
                else:
                    ls = line.split("\"")
                    arm_len[cc] = [float(ls[3])]
        elif "rgb(102,102,102)" in line:
            if "linearGradient" in line:
                pass
            elif "offset=" in line:
                pass
            else:
                if cc in arm_len:
                    ls = line.split("\"")
                    arm_len[cc].append(float(ls[3]))
                else:
                    ls = line.split("\"")
                    arm_len[cc] = [float(ls[3])]
        else:
            continue

    ilog.vlprint(f"arm_len dict {arm_len}", 8, logging=False)
    ilog.vlprint(f"arm_len.values() = {arm_len.values()}", 8, logging=False)
    chrom_scale = stat.mean([sum(sublist) for sublist in arm_len.values()])
    spec_res.chrom_scale = chrom_scale
    spec_res.calc_int_height()

    return templ, block_pos, arm_len


def gen_chrom_inst(templ, block_pos, arm_len, species_res, args):
    ilog.vlprint("Generating chromosome instances", 3)
    len_bp_real = species_res.len_bp_real
    len_bp_real = {int(k):v for k,v in len_bp_real.items()}
    species = species_res.spec_name
    int_height = species_res.int_height
    txt_displ = species_res.txt_displ
    corr_val = species_res.corr_val
    single_fig = species_res.single_fig
    ratios = []
    chrom_ins = []
    one_arm = [
        "Saccharomyces cerevisiae",
        "Caenorhabditis elegans",
        "Arabidopsis thaliana",
        "Danio rerio",
        "Drosophila melanogaster",
        "Escherichia coli"
        ]

    ilog.vlprint(f"block_pos", 6, logging=False)
    for _chr in block_pos:
        if species in one_arm :
            chr_type = "one_arm"
            ilog.vlprint(f"Chromosome {_chr}", 6, logging=False)
            ilog.vlprint(f"Arm-length {arm_len[_chr][0]}", 6, logging=False)
            ilog.vlprint(f"len_bp_real {len_bp_real[_chr]}", 6, logging=False)
            arm_length = arm_len[_chr][0] - corr_val
            cr_in = Chromosome(
                                _chr,
                                block_pos[_chr],
                                arm_length,
                                None,
                                len_bp_real[_chr],
                                templ,
                                chr_type,
                                int_height,
                                txt_displ,
                                single_fig,
                                args
                                )
        else:
            try:
                arm_two = arm_len[_chr][1]
            except IndexError:
                arm_two = 0
            chr_type = "two_arms"
            cr_in = Chromosome(
                                _chr,
                                block_pos[_chr],
                                arm_len[_chr][0],
                                arm_two,
                                len_bp_real[_chr],
                                templ,
                                chr_type,
                                int_height,
                                txt_displ,
                                single_fig,
                                args
                                )
        chrom_ins.append(cr_in)
        ratio = len_bp_real[_chr]/sum(arm_len[_chr])
        ratios.append(ratio)
    return chrom_ins


def gen_single_preambl(single_fig_dct, templ, prev, pas):
    ind_preamb = []
    transl = {"templ": templ, "prev": prev.block}
    for slice in single_fig_dct["ind_preamb"]:
        if slice[0] == "templ":
            x = pas
        else:
            x = 0
        if len(slice) == 3:
            if all(slice):
                ind_preamb.extend(transl[slice[0]][x+slice[1]: x+slice[2]])
            elif slice[1] == False:
                ind_preamb.extend(transl[slice[0]][: x+slice[2]])
            elif slice[2] == False:
                ind_preamb.extend(transl[slice[0]][x+slice[1]:])
            else:
                ilog.vlprint(
                    "WARNING: Unexpected option in gen_single_preambl", 0
                    )
        else:
            ind_preamb.append(transl[slice[0]][x+slice[1]])
    return ind_preamb


def draw_individual_chroms(
    preamble, chrom_ins, templ, templ_mod, css_styles, species, outdir, args, prefix=""
    ):
    "Draw ideograms of individual chromosomes containing integrations"
    ilog.vlprint(
        "Drawing ideograms of individual chromosomes containing integrations", 1
    )

    prev = None
    pas = preamble[species]

    for chr_obj in chrom_ins:
        if chr_obj.name == len(chrom_ins):
            chr_obj.last = True
        else:
            chr_obj.last = False

    os.mkdir(os.path.join(outdir, f"{prefix}single_chrom"))
    single_chroms = {}
    for chr_obj in chrom_ins:
        if chr_obj.last:
            co = chr_obj.single_fig["last_cutoff"]
            ind_preamb = gen_single_preambl(sfd, templ, prev, pas)
            for lin_idx in range(len(chr_obj.block) - co):
                ind_preamb.append(chr_obj.block[lin_idx])
            if co > 0:
                ind_preamb.append("</svg>")
            single_chroms[chr_obj.name] = ind_preamb
        elif prev:
            sfd = chr_obj.single_fig
            end = sfd["end"]
            ind_preamb = gen_single_preambl(sfd, templ, prev, pas)
            for lin_idx in range(len(chr_obj.block) - 8):
                ind_preamb.append(chr_obj.block[lin_idx])
            ind_preamb.append(end)
            single_chroms[chr_obj.name] = ind_preamb
        else:
            first_block = []
            end = chr_obj.single_fig["first_end"]
            for li in range(preamble[species] + len(css_styles)):
                templ_mod[li] = templ_mod[li].replace("opacity: 0", "opacity: 1")
                first_block.append(templ_mod[li])
            for lin_idx in range(len(chr_obj.block) - 8):
                first_block.append(chr_obj.block[lin_idx])
            first_block.append(end)
            single_chroms[chr_obj.name] = first_block
        prev = chr_obj

    ins_by_name = {obj.name:obj for obj in chrom_ins}
    for chrom in single_chroms:
        mcf = float(ins_by_name[chrom].mcf)
        fs = float(ins_by_name[chrom].font_size)/2
        prev_y = 0
        lime_line = ""
        fname = str(chrom) + ".svg"
        fpath = os.path.join(outdir, f"{prefix}single_chrom", fname)
        new_txt = {}
        # y = int()
        labl_cor = end = chr_obj.single_fig["labl_cor"]
        with open(fpath, "w") as ideo:
            for line in single_chroms[chrom]:
                if args.color_ints in line:
                    y1 = float(line.split("y=\"")[1].split("\"")[0])
                    y2 = y1
                    lime_line = f"                <line stroke=\"black\" stroke-width=\"{mcf/3}\" class=\"Integration\" x1=\"6\" x2=\"24\" y1=\"{y1}\" y2=\"{y2}\"/>\n"
                elif "blue" in line:
                    y1 = float(line.split("y=\"")[1].split("\"")[0])
                    y2 = y1
                    lime_line = f"                <line stroke=\"black\" stroke-width=\"{mcf/3}\" class=\"Integration\" x1=\"6\" x2=\"24\" y1=\"{y1}\" y2=\"{y2}\"/>\n"
                elif "class=\"int_txt\"" in line:
                    if "Position" in line:
                        new_txt["Pos"] = line.split(">Position [bp]: ")[1].split("</")[0]
                    elif "Features:" in line:
                        new_txt["Feat"] = line.split("Features: ")[1].split("</")[0]
                    else:
                        pass
                    if len(new_txt) == 2:
                        # mcf = float(ins_by_name[chrom].mcf)
                        # fs = float(ins_by_name[chrom].font_size)/2
                        dist_fit = y2/mcf - prev_y/mcf
                        if dist_fit <= fs:
                            y2_old = y2
                            y2 = prev_y + fs * mcf # fs - dist_fit + 
                            lime_line = lime_line.replace(
                                f"y2=\"{y2_old}\"", f"y2=\"{y2}\""
                                )
                            # adapt stroke-width according to angle of line
                            ankath, gegkath = y2-y1, 18 # x2-x1
                            hypot = math.sqrt(ankath**2 + gegkath**2)
                            sin_a = gegkath / hypot
                            angle_a = math.degrees(math.asin(sin_a))
                            cor_sw = (mcf/2) * angle_a/90
                            # cor_sw = (y2-y1)/mcf
                            lime_line = lime_line.replace(
                                f"stroke-width=\"{mcf/3}\"", f"stroke-width=\"{cor_sw}\""
                            )
                        int_txt = f'                <text class="int_txt" fill'\
                        f'="black" font-size="{fs}" transform="translate(0,0),'\
                        f'matrix(1 0 0 {mcf} 0 0)" x="32" y="'\
                        f'{(y2/mcf)+(fs/mcf)/2 + labl_cor}" >{new_txt["Pos"]}, '\
                        f'{new_txt["Feat"]}</text>\n'
                        ideo.write(lime_line)
                        ideo.write(int_txt)
                        new_txt = {}
                        prev_y = y2
                else:
                    ideo.write(line)

    int_chrom_ins = [str(ci.name) for ci in chrom_ins if len(ci.integrations) > 0]
    chrom_ideo_dir = os.scandir(os.path.join(outdir, f"{prefix}single_chrom"))
    for file in chrom_ideo_dir:
        if file.name.split(".")[0] not in int_chrom_ins:
            os.remove(os.path.abspath(file.path))


def ins_ints_write_out(templ, chrom_ins, outdir, outfile, line_spacer, chr_tr_dict, species, args, prefix=""):
    """Insert integrations and write to svg file"""
    ilog.vlprint("Inserting integrations and creating full ideogram", 1)
    ints_by_chr = defaultdict(list)
    ints_info = defaultdict(dict)

    features = dict()
    with open(os.path.join(outdir,f"{prefix}int_features.csv"), "r") as feat_csv:
            lc = 0
            for row in feat_csv:
                lc += 1
                if lc == 1:
                    pass
                else:
                    ls = row.strip().split(",")
                    if len(ls) > 2:
                        ls = [ls[0], ";".join(ls[1:])]
                    features[ls[0]] = ls[1]

    for intg in features: 
            feature = features[intg]
            col = feature.split("\t")
            req_info = [0, 1, 13, 14]
            col = [col[i] for i in range(len(col)) if i in req_info]
            if len(col) > 1:
                col[1] = f"({col[1]}) "
            if len(col) > 3:
                col[3] = f"{col[3]}"
            feature = " ".join(col)
            features[intg] = feature

    with open(os.path.join(outdir,f"{prefix}Integration_Report.csv"), "r") as int_csv:
        for line in int_csv:
            if line.startswith("ID"):
                pass
            else:
                ls = line.split(",")
                try:
                    chr_key = ls[2].strip()
                    chr_prefixes = ["NC", "NT", "CM", "BK", "BX"]
                    if chr_key not in chr_tr_dict:
                        ilog.vlprint(f"{chr_key} not in chr_tr_dict", 6)
                        if chr_key[:2] in chr_prefixes:
                            chr_key = chr_key.split(".")[0]
                            ilog.vlprint(f"Trying with {chr_key}", 6)
                    int_chr, int_pos_bp = chr_tr_dict[chr_key], int(ls[3])
                    int_id, int_name = ls[0], ls[1]
                    ints_by_chr[int_chr].append(int_pos_bp)
                    ints_info[int_chr][int_pos_bp] = {
                                                "int_id": int_id,
                                                "int_name": int_name,
                                                "contig_name": ls[2].strip(),
                                                "chr_name": int_chr,
                                                "int_pos_bp": int_pos_bp,
                                                "loc_accuracy": ls[4],
                                                "read_support": ls[5],
                                                "mms": ls[6],
                                                "ts_cov": ls[7],
                                                "IMR": ls[8],
                                                "SSBR": ls[9],
                                                "features": features[int_name]
                                                }
                except KeyError as ke:
                    ilog.vlprint(ke, 1)
                    ilog.vlprint(
                        f"INFO: Contig {ls[2].strip()} can not be represented"
                        f" in species specific histogram. This does only affect"
                        f" figure creation, not reports.", 2
                        )
                    ilog.vlprint(
                        f"INFO: This message will typically appear for sequences"
                        f" that are not embedded in the major chromosomes (e.g."
                        f" mitochondria, alternative loci).", 2
                        )
                    ilog.vlprint(
                        f"INFO: If this message occurs with major chromosomes, "
                        f"please check if the species (--spec) and reference-as"
                        f"sembly you specified are matching. Most likely there "
                        f"is a mismatch", 2
                        )
                    
    
    chrom_svg_lens = []
    for int_lst in ints_by_chr:
        for chr_inst in chrom_ins:
            chrom_svg_lens.append(int(chr_inst.total_len_svg))
            if int_lst == chr_inst.name:
                chr_inst.integrations = ints_by_chr[int_lst]
                chr_inst.int_info = ints_info[chr_inst.name]
    max_chrom_len = 0 if len(chrom_svg_lens)==0 else int(max(chrom_svg_lens))
    
    for inst in chrom_ins:
        ilog.vlprint(f"{inst.name} {int(inst.bp_svg_ratio)}", 7, logging=False)
        ilog.vlprint(f"{len(inst.block)} {inst.len_bp_total}", 7, logging=False)
        ilog.vlprint(
            f"{inst.total_len_svg} => {int(inst.total_len_svg * inst.bp_svg_ratio)}",
            7, logging=False
            )
        ilog.vlprint(f"{inst.integrations}", 7, logging=False)
        inst.max_chrom_len = max_chrom_len + 5
        if inst.bottom_aln:
            inst.transl = inst.transl - inst.max_chrom_len
        inst.insert_int(line_spacer)
    
    preamble = {
        "Saccharomyces cerevisiae": 58,
        "Homo sapiens": 139,
        "Mus musculus": 112,
        "Rattus norvegicus": 103,
        "Arabidopsis thaliana": 58,
        "Caenorhabditis elegans": 58,
        "Danio rerio": 58,
        "Drosophila melanogaster": 49,
        "Escherichia coli": int()
        }

    css_styles = [
                    "  <style>\n",
                    "    .int_txt {\n",
                    "      display: inline;\n",
                    "      text-align: center;\n",
                    "      color: none;\n",
                    "      opacity: 0\n",
                    "    }\n",
                    "    .Integration:hover {\n",
                    "      fill: red;\n",
                    "      stroke-color: red;"
                    "      color: red;\n",
                    "      pointer-events: visibleFill;\n",
                    "    }\n",
                    "    .Integration:hover ~ .int_txt {\n",
                    "      fill: black;\n",
                    "      color: black;\n",
                    "      pointer-events: visibleFill;\n",
                    "      opacity: 1\n",
                    "    }\n", 
                    "  </style>\n",
                ]

    lc = 0
    idx = 3
    for li in templ:
        lc += 1
        if "<desc>Ideogram SVG</desc>" in li:
            idx = lc
    templ_mod = [*templ[:idx], *css_styles[:], *templ[idx:]]
    
    # draw full ideogram
    with open(os.path.join(outdir, outfile), "w") as ideo:
        for li in range(preamble[species] + len(css_styles)):
            ideo.write(templ_mod[li])
        for chr_obj in chrom_ins:
            for lin_idx in range(len(chr_obj.block)):
                ideo.write(chr_obj.block[lin_idx])
    
    draw_individual_chroms(
        preamble, chrom_ins, templ, templ_mod, css_styles, species, outdir, args, prefix=prefix
        )


def gen_convertable(gradient_svg, species):
    """Generate a modified svg without gradients that can be converted to png"""
    ilog.vlprint("Generating convertable svg", 3)
    color_gradients = dict()
    tf_ct = 0
    transl_ct = 0
    transl = 0
    heigth = 0
    width = 12
    cur_transl = 0
    prev_transl = 0
    total_transl = 0
    svg_cycle = False
    block = False
    g_ct = 0
    corr_dict = {
        "Homo sapiens": 0,
        "Mus musculus": 0,
        "Saccharomyces cerevisiae": 24,
        "Arabidopsis thaliana": 35,
        "Caenorhabditis elegans": 35,
        "Rattus norvegicus": 0,
        "Drosophila melanogaster": 0,
        "Danio rerio": 24,
        "Escherichia coli": int()
        }
    
    convertable = open("ideo4png.svg", "w")
    with open(gradient_svg, "r") as gsvg:
        for line in gsvg:
            if line.startswith("<svg"):
                height = line.split("=")[1].split(" ")[0]
            elif "<linearGradient id=" in line:
                ls = line.strip().split("=")
                grad_name = ls[-1][1:-2]
                rgb_col = next(gsvg).strip().split("=")[2].split(" ")[0]
                color_gradients[grad_name] = rgb_col
            elif "<rect fill=\"url" in line:
                block = True
                ls = line.split("=")
                color_str, height_str = ls[1].split(" ")
                color_kw = color_str[6:-2]
                if color_kw in color_gradients:
                    ls.pop(1)
                    new_ins = " ".join([color_gradients[color_kw], height_str])
                    ls.insert(1, new_ins)
                    line = "=".join(ls)
            elif "<g transform=\"translate" in line:
                if "<g transform=\"translate(0," in line:
                    pass
                else:
                    tf_ct += 1
                    if tf_ct > 1:
                        prev_transl = cur_transl
                        transl_ct += 1
                        cur_transl = int(line.split("translate(")[1].split(",")[0])
                        if transl_ct == 2:
                            transl += cur_transl - corr_dict[species]
                        else:
                            transl += cur_transl - prev_transl
                        indent = line.split("<")[0]
                        line = f"{indent}<svg x=\"{transl}\" y=\"0\" width=\"12\" height={height}>\n"
            elif  "class=\"int_txt" in line:
                line = ""
            elif  "class=\"Integration" in line:
                g_ct -= 1
            elif "</g>" in line:
                if block:
                    g_ct += 1
                    if g_ct == 4:
                        indent = line.split("<")[0]
                        line = f"{indent}</svg>\n"
                        g_ct = 0
                        block = False
                else:
                    pass
            else:
                 pass
            convertable.write(line)
    convertable.close()


def run_spec_ideo(res_dir, outdir, outfile, line_spacer, species, args, prefix=""):
    """Creation of species specific ideogram with integration sites"""
    ilog.verbosity = args.verbosity
    report_time("generation of species specific ideograms")
    ilog.vlprint("Running species specific ideogram generation functions.", 3)
    spec_res = generate_spec_resource(species, res_dir)
    svg_src = os.path.join(res_dir, spec_res.svg_src)
    pss = parse_svg_src(svg_src, spec_res)
    gci = gen_chrom_inst(
        pss[0],
        pss[1],
        pss[2],
        spec_res,
        args
        )
    ins_ints_write_out(
        pss[0],
        gci,
        outdir,
        outfile,
        line_spacer,
        spec_res.chr_tr_dict,
        species,
        args,
        prefix=prefix,
        )
    gen_convertable(outfile, species)
    report_time("Generation of species specific ideograms", end=True)


if __name__ == "__main__":
    import sys
    svg_src = sys.argv[1]
    outdir = sys.argv[2]
    outfile = sys.argv[3]
    line_spacer = sys.argv[4]
    species = sys.argv[5]

    run_spec_ideo(svg_src, outdir, outfile, line_spacer, species, sys.argv)
