# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import csv
import os

from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import ImageReader
from reportlab.lib.units import cm
from reportlab.lib.colors import *
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM, renderPDF
from reportlab.platypus import *
from reportlab.lib.pagesizes import A4 # letter, A4 is default, landscape,
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER, TA_JUSTIFY
from reportlab.pdfbase.pdfmetrics import stringWidth

from il_logging import Ilogger


ilog = Ilogger()
ilog.module = __name__

def write_pdf(_Nx, species=None, avaclust=False):
    """Creates a pdf report"""
    ilog.vlprint("Creating pdf report", 1)
    fig_ct = 0
    tab_ct = 0
    # Prepare PDF
    styles = getSampleStyleSheet()
    # print(styles.__dict__)
    pdf = SimpleDocTemplate("Integration_summary_report.pdf", pagesize=A4)
    elements = []

    # # Generate Headline
    # prep logo
    res_dir = os.path.join(os.path.split(os.path.abspath(__file__))[0], "il_resources")
    logo = Image(os.path.join(res_dir, "il_logo.png")) # , 1.5*cm, 26.5*cm
    logo.drawWidth, logo.drawHeight = 1.8*cm, 1*cm
    # prep Title
    title_style = styles["Heading1"]
    # print(title_style.__dict__)
    title_style.fontSize = 48
    main_Title = Paragraph("Summary report", title_style)
    # set headline
    header_table = Table([[logo, main_Title]], [2*cm, 15*cm], spaceAfter=54)
    header_table.wrap(16 * cm, 15 * cm)
    elements.append(header_table)

    # # Generate Run-Info table
    # make T1 caption
    tab_ct += 1
    title_style.fontSize = 18
    title_table1 = Paragraph(f"Table {tab_ct}: Run Information", title_style)
    elements.append(title_table1)
    # make Table1
    with open('run_info.csv', "r") as csvfile:
        data = list(csv.reader(csvfile))
        for i in range(len(data)):
            if data[i][0] == "Output directory":
                pass
            elif data[i][0] == "sys_version":
                sys_version = data[i][1]
            else:
                data[i][1] = os.path.split(data[i][1])[-1]
        for i in range(len(data)):
            if data[i][0] == "OS":
                data[i][1] = data[i][1] + "_" + sys_version

    filt = ["sys_release", "sys_version", "num_CPUs", "machine", "processor", ]
    data = [entry for entry in data if entry[0] not in filt]
    for col_lst in range(len(data)):
        for entry in range(len(data[col_lst])):
            data[col_lst][entry] = Paragraph(data[col_lst][entry], styles['Normal'])

    inf_table = Table(data, [6*cm, 10*cm], spaceAfter=36)
    inf_tab_style = TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LINEBELOW', (0, 0), (-1, 0), 1, black),
        ('ALIGN', (0, 0), (-1, 0), 'LEFT'), # CENTER
        ('FONT', (0, 0), (-1, 0), "Helvetica-Bold", 13),
        ('GRID', (0, 0), (-1, -1), 1, black),
        ('ALIGN', (3, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (2, 1), (2, -1), 'RIGHT'),
        ])
    inf_table.setStyle(inf_tab_style)
    elements.append(inf_table)

    def resize_image(img_path, width=15 * cm):
        img = ImageReader(img_path)
        orig_width, orig_height = img.getSize()
        aspect_ratio = float(orig_height) / float(orig_width)
        return Image(img_path, width=width, height=(width * aspect_ratio))

    # # Plot Relative positions of Integration locations Figure
    fig_ct += 1
    title_fig1 = Paragraph(f"Fig. {fig_ct}: Relative positions of Integrations in contigs", title_style)
    # elements.append(title_fig1)
    svg_obj = svg2rlg("relative_integration_locations.svg")
    renderPM.drawToFile(svg_obj, "relative_integration_locations.png", fmt="PNG") # , dpi=72
    fig1_plot = resize_image("relative_integration_locations.png")# 2.5*cm, 1.5*cm, width=15*cm,
    fig1_plot.spaceAfter = 36
    f1_compound = KeepTogether([title_fig1, fig1_plot])
    # canvas.drawString(2*cm, 5.5*cm, "Figure 1: Relative positions of Integrations in Contigs")
    elements.append(f1_compound)

    # # Generate Integrations Table
    # make T2 caption + subtext
    tab_ct += 1
    title_table2 = Paragraph(f"Table {tab_ct}: Integration sites (read to reference alignm.) ", title_style)
    subtext_style = styles["Normal"]
    subtext_style.fontSize = 12
    subtext_style.spaceAfter = 12
    subtext_style.alignment = TA_CENTER
    subtext_t2 = Paragraph("(for names of supporting reads see \"Integration_Report.csv\")", subtext_style)
    elements.append(title_table2)
    elements.append(subtext_t2)

    with open('Integration_Report.csv', "r") as csvfile:
        data = list(csv.reader(csvfile))
    d_short = data.copy()
    for row in d_short:
        row.pop()

    int_tab_style = TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LINEBELOW', (0, 0), (-1, 0), 1, black),
        ('ALIGN', (0, 0), (-4, 0), 'LEFT'),
        ('ALIGN', (-4, 0), (-1, 0), 'CENTER'),
        ('FONT', (0, 0), (-1, 0), "Helvetica-Bold", 10),
        ('BACKGROUND', (0, 0), (-1, 0), lightblue),
        ('GRID', (0, 0), (-1, -1), 1, black),
        ('ALIGN', (3, 1), (-1, -1), 'CENTER'),
        ('ALIGN', (2, 1), (2, -1), 'RIGHT'),
        ])

    # dynamically set column-widths, bounded by min/max values
    def get_cw(col, _data, font="Helvetica", font_size=12, min_width={}):
        """Calculate required width in pixels for a column"""
        col_items = []
        for row_i in range(len(_data)):
            col_items.append(stringWidth(str(_data[row_i][col]), font, font_size))
        try:
            if max(col_items) > min_width[col]:
                return max(col_items)
            else:
                return min_width[col]
        except KeyError:
            print("KeyError in get_cw")
            return max(col_items)

    min_w_t2 = {0: 25, 1: 60, 2: 60, 3: 30, 4: 40, 5: 30, 6: 30, 7: 30, 8: 30}
    max_w_t2 = {0: 60, 1: 120, 2: 120, 3: 50, 4: 100, 5: 60, 6: 60, 7: 60, 8: 60}
    c_widths = [get_cw(col, d_short, min_width=min_w_t2) for col in min_w_t2]
    if sum(c_widths) > sum(max_w_t2.values()):
        c_widths = list(max_w_t2.values())
    col_widths = [
        c_widths[0],
        c_widths[1],
        c_widths[2],
        c_widths[3],
        c_widths[4], 
        c_widths[5],
        c_widths[6],
        c_widths[7],
        c_widths[8],
        ]

    int_table = Table(d_short, colWidths=col_widths, spaceAfter=36, repeatRows=1)
    int_table.setStyle(int_tab_style)
    elements.append(int_table)

    # # Plot Ints/Chromosome Figure
    fig_ct += 1
    title_fig2 = Paragraph(f"Fig. {fig_ct}: Integrations/contig (int.-positive contigs)", title_style)
    fig2_plot = resize_image("ints_per_chrom.png")# 2.5*cm, 1.5*cm, width=15*cm,
    fig2_plot.spaceAfter = 36
    f2_compound = KeepTogether([title_fig2, fig2_plot])
    elements.append(f2_compound)

    # # Plot Ints/Chromosome Figure
    fig_ct += 1
    title_fig3 = Paragraph(f"Fig. {fig_ct}: Integrations/contig (contigs N{_Nx})", title_style)
    fig3_plot = resize_image("ints_per_chrom_incl0.png")# 2.5*cm, 1.5*cm, width=15*cm,
    fig3_plot.spaceAfter = 36
    f3_compound = KeepTogether([title_fig3, fig3_plot])
    elements.append(f3_compound)

    if avaclust:
        try:
            # Plot results from all-vs-all read clustering
            # Create Table listing clusters and reads
            tab_ct += 1
            title_table3 = Paragraph(
                f"Table {tab_ct}: Clusters of reads with high identity "
                f"(all-vs-all alignm.)", title_style
                )
            elements.append(title_table3)

            with open('ava_read_cluster.csv', "r") as csvfile:
                data_t3 = list(csv.reader(csvfile))

            for row_t3 in range(len(data_t3)):
                if row_t3 == 0:
                    pass
                else:
                    for col_t3 in range(len(data_t3[row_t3])):
                        style = styles['Normal']
                        if col_t3 == 2:
                            data_t3[row_t3][col_t3] = data_t3[row_t3][col_t3].replace(";", " // ")
                            style.fontSize = 6
                        else:
                            style.fontSize = 12
                        data_t3[row_t3][col_t3] = Paragraph(data_t3[row_t3][col_t3], styles['Normal'])

            clust_tab_style = TableStyle([
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('LINEBELOW', (0, 0), (-1, 0), 1, black),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONT', (0, 0), (-1, 0), "Helvetica-Bold", 13),
                ('BACKGROUND', (0, 0), (-1, 0), red),
                ('GRID', (0, 0), (-1, -1), 1, black),
                ])

            min_w_t3 = {0: 60, 1: 60, 2: 330, } # 3: 65, 4: 55,
            max_w_t3 = {0: 60, 1: 60, 2: 330, } # 3: 65, 4: 55,
            c_widths_t3 = [get_cw(col, data_t3, min_width=min_w_t3) for col in min_w_t3]
            if sum(c_widths_t3) > sum(max_w_t3.values()):
                c_widths_t3 = list(max_w_t3.values())
            col_widths_t3 = [c_widths_t3[0], c_widths_t3[1], c_widths_t3[2], ]

            table3 = Table(data_t3, colWidths=col_widths_t3, spaceAfter=36, repeatRows=1)
            table3.setStyle(clust_tab_style)
            elements.append(table3)


            # # # Integrate Cluster Heatmaps
            # # il-ava-clust
            # make Fig4a caption + subtext
            fig_ct += 1
            title_fig4a = Paragraph(f"Fig. {fig_ct}a: Read groups from all-vs-all alignment", title_style)
            subtext_style.spaceAfter = 12
            subtext_style.fontSize = 10
            subtext_f4a = Paragraph("clustered by bitscore with il_ava_clust; qcovus values shown", subtext_style)

            fig4a_plot = resize_image("read_cluster_il-ava-clust.png")
            fig4a_plot.spaceAfter = 36
            f4a_compound = KeepTogether([title_fig4a, subtext_f4a, fig4a_plot])
            elements.append(f4a_compound)

            # # scipy-ava-clust + dendrogram
            # make Fig4b caption + subtext
            title_fig4b = Paragraph(f"Fig. {fig_ct}b: Read groups from all-vs-all alignment", title_style)
            subtext_style.spaceAfter = False
            subtext_style.fontSize = 10
            subtext_f4b = Paragraph("clustered by qcovus with scipy; qcovus values shown;"
                                    "bs_cluster < min_read_supp excluded;", subtext_style)
            subtext_style.spaceAfter = 12

            fig4b_plot = resize_image("read_cluster_ava_scipy.png", width=17 * cm)
            fig4b_plot.spaceAfter = 36
            f4b_compound = KeepTogether([title_fig4b, subtext_f4b, fig4b_plot])
            elements.append(f4b_compound)
        except:
            title_table3 = Paragraph(f"Read Clustering", title_style)
            elements.append(title_table3)
            subtext_style.spaceAfter = False
            subtext_style.fontSize = 10
            subtext_f4b = Paragraph("Read clustering could not be performed", subtext_style)

    if species:
        try:
            fig_ct += 1
            spec_ideo_svg = svg2rlg("ideo4png.svg")
            renderPM.drawToFile(spec_ideo_svg, "ideogram.png", fmt="PNG")
            title_fig5 = Paragraph(
                f"Fig. {fig_ct}: Species specific ideogram", title_style
                )
            ideo_width = {
            "Drosophila melanogaster": 10,
            "Homo sapiens": 17,
            "Saccharomyces cerevisiae": 19,
            "Danio rerio": 19,
            "Mus musculus": 19,
            "Escherichia coli": 19,
            "Arabidopsis thaliana": 19,
            "Rattus norvegicus": 19,
            "Caenorhabditis elegans": 19,
            }
            fig5_plot = resize_image(
                "ideogram.png", width=ideo_width[species] * cm
                )
            fig5_plot.spaceAfter = 36
            f5_compound = KeepTogether([title_fig5, fig5_plot])
            elements.append(f5_compound)
        except Exception as e:
            ilog.vlprint(f"WARNING: Error in il_pdf species ideo {e}", 0)
            ilog.vlprint(
                f"WARNING: Species specific ideogram could not be included in"
                f" pdf report, likely due to missing ideogram svg file.", 0
                )

        # Plot features if available
        # Create Table listing clusters and reads
        tab_ct += 1
        title_table4 = Paragraph(f"Table {tab_ct}: Genomic features at Integration site", title_style)
        elements.append(title_table4)

        with open('int_features.csv', "r") as csvfile:
            lc = 0
            data_t4 = dict()
            for row_t4 in csvfile:
                lc += 1
                if row_t4 == 1:
                    pass
                else:
                    ls = row_t4.strip().split(",")
                    if len(ls) > 2:
                        ls = [ls[0], ";".join(ls[1:])]
                    data_t4[ls[0]] = ls[1]

        for intg in data_t4: 
            feature = data_t4[intg]
            style = styles['Normal']
            col = feature.split("\t")
            req_info = [0, 1, 13, 14]
            col = [col[i] for i in range(len(col)) if i in req_info]
            if len(col) > 1:
                col[1] = f"({col[1]}): "
            if len(col) > 3:
                col[3] = f"{col[3]};"
            feature = " ".join(col)
            style.fontSize = 11
            data_t4[intg] = Paragraph(feature, styles['Normal'])
        data_t4 = list(data_t4.items())

        feat_tab_style = TableStyle([
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LINEBELOW', (0, 0), (-1, 0), 1, black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONT', (0, 0), (-1, 0), "Helvetica-Bold", 13),
            ('BACKGROUND', (0, 0), (-1, 0), orange),
            ('GRID', (0, 0), (-1, -1), 1, black),
        ])

        def get_cw(col, _data, font="Helvetica", font_size=12, min_width={}):
            """Calculate required width in pixels for a column"""
            col_items = []
            for row_i in range(len(_data)):
                col_items.append(stringWidth(str(_data[row_i][col]), font, font_size))
            try:
                if max(col_items) > min_width[col]:
                    return max(col_items)
                else:
                    return min_width[col]
            except KeyError:
                return max(col_items)

        min_w_t4 = {0: 120, 1: 120}
        max_w_t4 = {0: 180, 1: 275}
        try:
            c_widths_t4 = [get_cw(col, data_t4, min_width=min_w_t4) for col in min_w_t4]
        except IndexError as ie:
            print("Dynamic assignment of column widths failed for feature table")
            c_widths_t4 = list(max_w_t4.values())
        if sum(c_widths_t4) > sum(max_w_t4.values()):
            c_widths_t4 = list(max_w_t4.values())
        col_widths_t4 = [c_widths_t4[0], c_widths_t4[1]]

        table4 = Table(data_t4, colWidths=col_widths_t4, spaceAfter=36, repeatRows=1)
        table4.setStyle(feat_tab_style)
        elements.append(table4)

    pdf.build(elements)


if __name__ == "__main__":
    Nx = 98
    species = "human" # False
    write_pdf(Nx, species=species)
