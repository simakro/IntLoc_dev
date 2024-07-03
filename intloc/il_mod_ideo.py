# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import gzip
import os
import re
from collections import Counter, defaultdict

from il_logging import Ilogger
from il_spec_ideo import SpeciesResources, generate_spec_resource


ilog = Ilogger()
ilog.module = __name__


def gene_density(features, spec_res, tfeat="gene", tclass="protein_coding"):
    features = gzip.open(features, "r")
    bin_size = spec_res.gene_dens["bin_size"]
    feat_class_by_chrom = defaultdict(list)
    feat_class_ct = 0
    feat_count = Counter()
    gene_count = Counter()
    cols = [
            "feature",
            "class",
            "assembly",
            "assembly_unit",
            "seq_type",
            "chromosome",
            "genomic_accession",
            "start",
            "end",
            "strand",
            "product_accession",
            "non-redundant_refseq",
            "related_accession",
            "name",
            "symbol",
            "GeneID",
            "locus_tag",
            "feature_interval_length",
            "product_length",
            "attributes"
            ]
    all_classes = set()
    all_features = set()
    for line in features:
        line = line.decode("utf-8")
        if line.startswith("#"):
            pass
        else:
            ls = line.split("\t")
            feat_dct = dict(zip(cols, ls))
            feat_count[ls[0]] += 1
            all_features.add(feat_dct["feature"])
            all_classes.add(feat_dct["class"])
            if feat_dct["feature"] == "gene":
                gene_count[feat_dct["class"]] += 1
            if feat_dct["feature"] == tfeat and feat_dct["class"] == tclass:
                feat_class_by_chrom[feat_dct["chromosome"]].append(feat_dct)
                feat_class_ct += 1

    ilog.vlprint(f"Available features {all_features} for feature-ideo", 5)
    ilog.vlprint(f"Available classes {all_classes} for feature-ideo", 5)

    feat_class_by_chrom = {
        str(spec_res.chr_tr_dict[str(k)]):v for k,v in feat_class_by_chrom.items() if len(k) > 0
        }
    for chrom in spec_res.len_bp_real:
        if chrom not in feat_class_by_chrom:
            feat_class_by_chrom[chrom] = []
    
    chroms_bins = {}
    for chrom in feat_class_by_chrom:
        try:
            chr = str(spec_res.chr_tr_dict[chrom])
            chr_len = spec_res.len_bp_real[chr]
            num_bins = range(int(chr_len / bin_size)+1)
            chr_bins = {((bin*bin_size-(bin_size-1)),bin*bin_size):[] for bin in num_bins}
            chr_bins.pop((-(bin_size-1), 0))
            # adjust bins for chr-overhang not included in last bin, by either
            # extending last bin or adding another bin depending on overhangsize
            last_bin = list(chr_bins.keys())[-1]
            chr_overhang = chr_len - last_bin[1]
            if chr_overhang < 0.5*bin_size:
                chr_bins.pop(last_bin)
                mod_lbin = (last_bin[0], last_bin[1]+chr_overhang)
                chr_bins[mod_lbin] = []
            else:
                new_bin = (last_bin[1]+1, last_bin[1]+chr_overhang)
                chr_bins[new_bin] = []
            chr_bins = dict(sorted(chr_bins.items()))
            feat_lst = sorted(feat_class_by_chrom[chrom], key=lambda x: x["start"])
            for feat in feat_lst:
                for fbin in chr_bins:
                    if fbin[0] <= int(feat["start"]) <= fbin[1]:
                        chr_bins[fbin].append(feat)
                        if int(feat["end"]) <= fbin[1]:
                            break
                    elif fbin[0] <= int(feat["end"]) <= fbin[1]:
                        chr_bins[fbin].append(feat)
                        break
                    else:
                        pass
            chroms_bins[chr] = chr_bins
        except ValueError as ve:
            ilog.vlprint(f"Gene-density for chrom {ve} can not be plotted", 1)
        except KeyError as ke:
            ilog.vlprint(
                f"Chrom {ke} could not be found in the reference. This is"
                f" usually related to the mitochondrion and no reason to worry."
                ,1
                )

    chrom_gene_dens = {}
    for chrom in chroms_bins:
        gd_dist = []
        chr_len = spec_res.len_bp_real[chr]
        for fbin in chroms_bins[chrom]:
            feat_symbols = [feat["symbol"] for feat in chroms_bins[chrom][fbin]]
            feat_set = set(feat_symbols)
            gd_dist.append(len(feat_set))
        chrom_gene_dens[chrom] = gd_dist
    return chrom_gene_dens


def spectral_gradient():
    """Set up a spectral color gradient"""
    rgb_curr = [0,0,0]
    rgb_dict = {}
    for n in range(1021):
        if n == 0:
            rgb_curr = [0,0,255]
            rgb_dict[n] = rgb_curr.copy()
        elif 0<n<256:
            rgb_curr[1] += 1
            rgb_dict[n] = rgb_curr.copy()
        elif 254<n<511:
            rgb_curr[2] -= 1
            rgb_dict[n] = rgb_curr.copy()
        elif 509<n<766:
            rgb_curr[0] += 1
            rgb_dict[n] = rgb_curr.copy()
        elif 764<n<1021:
            rgb_curr[1] -= 1
            rgb_dict[n] = rgb_curr.copy()
    return rgb_dict


def spectral_grad_mids(step_size=5):
    """Set up a spectral color gradient with emphasized mid-tones"""
    rgb_curr = [0,0,0]
    rgb_dict = {}
    for n in range(1021):
        if n == 0:
            rgb_curr = [0,0,255]
            rgb_dict[n] = rgb_curr.copy()
        elif 0<n<256:
            rgb_curr[1] += step_size
            rgb_dict[n] = rgb_curr.copy()
        elif 254<n<511:
            rgb_curr[2] -= 1
            rgb_dict[n] = rgb_curr.copy()
        elif 509<n<766:
            rgb_curr[0] += 1
            rgb_dict[n] = rgb_curr.copy()
        elif 764<n<1021:
            rgb_curr[1] -= step_size
            rgb_dict[n] = rgb_curr.copy()
    rgb_dict = {
        k:[val if val<256 else 255 for val in v] for k,v in rgb_dict.items()
        }
    return rgb_dict


def get_rgb_val(num, max_num, rgb_dict):
    min_incr = int(max_num)/(4*255)
    scaled_num = int(int(num)/min_incr)
    return rgb_dict[scaled_num]


def scan_single_fig_dir(outdir):
    ilog.vlprint(f"Scanning figure dir",3)
    sf_dir = os.path.join(outdir, "single_chrom")
    scan = os.scandir(sf_dir)
    sc_figs = [entry.path for entry in scan if os.path.split(entry)[1].split(".")[-1]=="svg"]
    return sc_figs


def scan_outdir(outdir):
    ilog.vlprint(f"Scanning outdir {outdir}",3)
    scan = os.scandir(outdir)
    sc_figs = [entry.path for entry in scan if os.path.split(entry)[1]=="ideogram.svg"]
    return sc_figs


def process_ideo_batch(
    sc_figs,
    rgb_dict,
    spec_res,
    gds_chroms,
    tfeat="gene",
    tclass="protein_coding",
    repl_label=True,
    single_chrom=True
    ):
    """Process a batch of template ideograms."""
    ilog.vlprint("Processing batch of template ideograms.", 3)
    max_gd_genome = max([num for lst_val in gds_chroms.values() for num in lst_val])
    for fig in sc_figs:
        if os.path.split(fig)[-1] == "ideogram.svg":
            single_chrom = False
        gf = write_grad_ideo(
            fig,
            gds_chroms,
            rgb_dict,
            spec_res,
            max_num_gen=max_gd_genome,
            single_chrom=single_chrom,
            tfeat=tfeat,
            tclass=tclass,
            )
        if repl_label:
            replace_marker_triangle(gf, spec_res)
    return max_gd_genome


def write_grad_ideo(
    svg_tmpl,
    gds_chroms,
    rgb_dict,
    spec_res,
    max_num_gen=False,
    single_chrom=True,
    tfeat="gene",
    tclass="protein_coding"
    ):
    """Write out ideograms with feature density information"""
    ilog.vlprint("Writing out ideograms with feature density information", 3)
    if max_num_gen:
        max_num = max_num_gen
    else:
        max_num = max(gd_values)
    with open(svg_tmpl, "r") as svg_in:
        template = svg_in.read().split("\n")
    idx = -1
    cp_fill = False
    offset = ""
    y = float()
    old_fill = []
    curr_chrom = ""
    path_lens = {}
    bin_size = float()
    for line in template:
        idx += 1
        if '<text font-family="Courier" font-size="12.000000px"' in line:
            curr_chrom = line.split(">")[1].split("<")[0]
            try:
                gd_values = gds_chroms[curr_chrom]
            except KeyError as ke1:
                ilog.vlprint(
                    f"No gene-density information was found for chromosome"
                    f" {ke1}. Attempting to retrieve info with chr-name synonym"
                    , 5
                    )
                try:
                    gd_values = gds_chroms[str(spec_res.chr_tr_dict[curr_chrom])]
                    ilog.vlprint(
                    f"Retrieval of gene-density information successful.", 5
                    )
                except KeyError as ke2:
                    ilog.vlprint(
                    f"No gene-density information was found for chromosome"
                    f" {ke1} or synonym {ke2}. Skipping chromosome."
                    , 5
                    )
        if 'class="Integration"' in line:
            template.pop(idx)
            template.insert(idx, re.sub('fill="[a-z]*"', 'fill="black"', line))
        if '<path id="' in line:
            ls = line.split("\"")
            key, val = ls[1], ls[-8].split(" ")[-4].split(",")[1]
            path_lens[key] = val
        if '<use xlink:href="url(' in line or '<use href=' in line:
            curr_path = line.split("(")[1].split(")")[0].split("#")[1]
            plen = float(path_lens[curr_path])
            bin_size = (plen)/len(gd_values)
        if "</g>" in line:
            if cp_fill:
                first_idx = old_fill[0]
                delta_len = len(gd_values)-len(old_fill)
                for _ in old_fill:
                    template.pop(first_idx) 
                for val in gd_values:
                    y += bin_size
                    rgb_val = get_rgb_val(val, max_num, rgb_dict)
                    nl = f'<rect fill="rgb{tuple(rgb_val)}" height="{bin_size}" width="12" stroke="none" x="0" y="{y}"/>'
                    template.insert(first_idx, f'{offset}{nl}')
                    first_idx += 1
                if delta_len < 0:
                    for _ in range(abs(delta_len)):
                        template.insert(first_idx, "\n")
                        first_idx += 1
                old_fill = []
            cp_fill = False
        elif cp_fill:
            offset = line.split("<")[0]
            old_fill.append(idx)
        elif "<g clip-path=" in line:
            y = 0 - bin_size
            cp_fill = True
        else:
            pass
    template = [line for line in template if len(line) > 1]
    if single_chrom:
        out_file = os.path.join(os.path.split(svg_tmpl)[0], f'{curr_chrom}_dens_{tclass}_{tfeat}.svg')
        with open(out_file, 'w') as out:
            out.write("\n".join(template))
    else:
        out_file = os.path.join(os.path.split(svg_tmpl)[0], f'ideogram_dens_{tclass}_{tfeat}.svg')
        with open(out_file, 'w') as out:
            out.write("\n".join(template))
    return out_file


def replace_marker_triangle(svg_tmpl, spec_res):
    """Change integration marker to triangle"""
    ilog.vlprint("Changing integration marker to triangle", 3)
    with open(svg_tmpl, "r") as svg_in:
        template = svg_in.read().split("\n")
    idx = -1
    scale = float(spec_res.txt_displ["mcf"])
    for line in template:
        idx += 1
        if 'class="Integration"' in line:
            try:
                x = int(line.split("x1=")[1].split(" ")[0][1:-1]) + 6
                y = float(line.split("y1=")[1].split(" ")[0][1:-1])
            except:
                x = 12
                y = float(line.split('y="')[1].split('"')[0])
            offset = line.split("<")[0]
            template.pop(idx)
            bc = [0, 0, 3, 3*scale, 3, -3*scale]
            pc = [x+bc[0], y+bc[1], x+bc[2], y+bc[3], x+bc[4], y+bc[5]]
            nl = f'{offset}<polygon class="int-site" fill="green" x="{x}" y="{y}" points="{pc[0]} {pc[1]}, {pc[2]} {pc[3]}, {pc[4]} {pc[5]}"/>'
            template.insert(idx, nl)
        if '<text class="int_txt"' in line:
            template.pop(idx)
            template.insert(idx, "")

    out_file = ".".join(svg_tmpl.split(".")[:-1]) + "_tr"  + ".svg"
    with open(out_file, 'w') as out:
        out.write("\n".join(template))


def write_grad_legend(rgb_dict, bin_size=1, label_interval=100, max_num=1020):
    """Prepare color gradient legend for feature density plots"""
    ilog.vlprint("Preparing color gradient legend for feature density plots", 3)
    preamble = (
"""<?xml version="1.0"?>
<svg height="1200" version="1.1" width="600" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
    <desc>Gradient</desc>
"""
    )
    labels = [n*label_interval for n in range(int(len(rgb_dict)/label_interval)+1)]
    values = range(len(rgb_dict))
    close_svg = ('</svg>')
    with open("Test.svg", "w") as out:
        out.write(preamble)
        y=10
        for val in values:
            y += 1
            rgb_val = rgb_dict[val]
            line = f'    <line stroke="rgb{tuple(rgb_val)}" stroke-width="{bin_size}" x1="100" x2="200" y1="{y}" y2="{y}"/>\n'
            out.write(line)
            if val in labels:
                min_incr = int(max_num)/(4*255)
                scaled_num = int(int(val)*min_incr)
                out.write(f'    <line stroke="black" stroke-width="{bin_size}" x1="200" x2="220" y1="{y}" y2="{y}"/>\n')
                out.write(f'    <text stroke="black" class="int_txt" fill="black" font-size="8.0" x="225" y="{y}" >{scaled_num}</text>\n')
        out.write(close_svg)


def run_mod_ideo(args, tfeat="gene", tclass="protein_coding"):
    """Run the mod_ideo module to derive customized ideograms from templates"""
    ilog.vlprint(f"Creating derived {tclass} {tfeat} density ideograms.",1)
    spec_res = generate_spec_resource(args.species, args.res_dir)
    _features = spec_res.feat_src
    rgb_dict = spectral_grad_mids(spec_res.gene_dens["grad_step_size"])
    gds_chroms = gene_density(
        _features, spec_res, tfeat=tfeat, tclass=tclass
        )
    sc_figs = scan_single_fig_dir(args.outdir)
    max_gd_genome = process_ideo_batch(
        sc_figs, rgb_dict, spec_res, gds_chroms, tfeat=tfeat, tclass=tclass
        )
    full_ideo = scan_outdir(args.outdir)
    process_ideo_batch(
        full_ideo, rgb_dict, spec_res, gds_chroms, single_chrom=False, tfeat=tfeat, tclass=tclass
        )
    return gds_chroms, rgb_dict, max_gd_genome


if __name__ == "__main__":
    import sys

    # def run_mod_ideo(species, res_dir, out_dir, _features, target_feat="", target_class=""):
    #     spec_res = generate_spec_resource(species, res_dir)
    #     rgb_dict = spectral_grad_mids(spec_res.gene_dens["grad_step_size"])
    #     gds_chroms = gene_density(_features, spec_res, target_feat=target_feat, target_class=target_class)
    #     sc_figs = scan_single_fig_dir(out_dir)
    #     max_gd_genome = process_ideo_batch(sc_figs, rgb_dict, spec_res, gds_chroms)
    #     full_ideo = scan_ideo_dir(out_dir)
    #     process_ideo_batch(full_ideo, rgb_dict, spec_res, gds_chroms, single_chrom=False)
    #     return gds_chroms, rgb_dict, max_gd_genome

    # _features = gzip.open(sys.argv[1], "r")
    # target_feat = sys.argv[2]
    # target_class = sys.argv[3]
    # species = sys.argv[4]
    # res_dir = "../intloc/il_resources"
    # out_dir = "/il_test/Test_results/intloc_out_2022-12-19_11-55-36_yeast_standard_reports"
    
    # gds, rgb, mgg = run_mod_ideo(
    #     species, res_dir, out_dir, _features, target_feat=target_feat, target_class=target_class
    #     )
    # write_grad_legend(rgb, bin_size=1, max_num=mgg)
