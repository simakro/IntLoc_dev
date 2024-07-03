# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import sys
import os
import random

svg = sys.argv[1]
genome_info = sys.argv[2]

with open(svg, "w") as svg:
    preamble = ("<?xml version=\"1.0\"?>\n"
    "<svg height=\"320\" version=\"1.1\" width=\"640\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n"
    "  <desc>Ideogram SVG</desc>\n"
    "  <defs>\n"
    "    <linearGradient id=\"arm\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(102,102,102)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(195,195,195)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(136,136,136)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(102,102,102)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\n"
    "    <linearGradient id=\"gneg\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(191,191,191)\" stop-opacity=\"1\"/>‚\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(255,255,255)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(255,255,255)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(191,191,191)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\n"
    "    <linearGradient id=\"gpos-100\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(0,0,0)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(127,127,127)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(0,0,0)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(0,0,0)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\n"
    "    <linearGradient id=\"gpos-33\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(127,127,127)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(212,212,212)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(170,170,170)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(127,127,127)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\n"
    "    <linearGradient id=\"gpos-66\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(76,76,76)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(178,178,178)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(102,102,102)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(76,76,76)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\n"
    "    <linearGradient id=\"gpos-75\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(51,51,51)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(161,161,161)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(68,68,68)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(51,51,51)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\n"
    "    <linearGradient id=\"acen\">\n"
    "      <stop offset=\"0\" stop-color=\"rgb(191,146,146)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.125\" stop-color=\"rgb(255,163,163)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"0.625\" stop-color=\"rgb(255,192,192)\" stop-opacity=\"1\"/>\n"
    "      <stop offset=\"1\" stop-color=\"rgb(191,146,146)\" stop-opacity=\"1\"/>\n"
    "    </linearGradient>\"\n"
    "  </defs>\n"
    "  <g transform=\"translate(0,0),matrix(1 0 0 1 0 0)\">\n"
    "    <g transform=\"translate(0,0),matrix(1 0 0 1 0 0)\">\n"
    "      <g transform=\"translate(20,20),matrix(1 0 0 1 0 16)\">\n"
    )

    svg.write(preamble)

    
    with open(genome_info, "r") as csv:
        cct = 0
        for line in csv:
            if line.startswith("chr_num"):
                pass
            else:
                cct += 1
                i = cct
                transl = cct * 35
                i = i+1
                clip_path_id = i*2
                path_id = clip_path_id -1
                ls = line.split(",")
                chrom = ls[3]
                p_arm = int(ls[1]) / 100000
                cen = int(ls[-2]) / 100000
                #cen_width = (int(ls[-1]) - int(ls [-2])) / 5000
                #q_arm = 2*p_arm + random.randint(10,200)
                chr_block = (
                f"        <clipPath id=\"ideo_{clip_path_id}\">\n"
                f"          <path id=\"ideo_{path_id}\" d=\"M0,5.05 C0,2.5 3,0 6,0 C9,0 12,2.5 12,5 L12,{p_arm} C12,{p_arm+ 2} 9,{p_arm + 5} 6 {p_arm + 5} C3,{p_arm + 5} 0,{p_arm + 2} 0,{p_arm} Z\" fill=\"rgb(240,240,255)\" stroke=\"rgb(0,0,0)\" stroke-width=\"1.000px\"/>\n" # M0,{p_arm + 10} C0,{p_arm + 7} 3,{p_arm + 5} 6,{p_arm + 5} C9,{p_arm + 5} 12,{p_arm + 7} 12,{p_arm + 10} L12,{q_arm} C12,{q_arm + 3} 9,{q_arm + 5} 6 {q_arm + 5} C3,{q_arm + 5} 0,{q_arm + 3} 0,{q_arm} Z\"
                f"        </clipPath>\n"  
                f"        <g transform=\"translate({transl},0),matrix(1 0 0 1 0 0)\">\n"
                f"          <g fill=\"rgb(32,32,32)\" fill-opacity=\"1\" transform=\"translate(3,-16),matrix(1 0 0 1 0 0)\">\n"
                f"            <text font-family=\"Courier\" font-size=\"12.000000px\" transform=\"matrix(1 0 0 1 0 8)\" x=\"0\" y=\"0\">{chrom}</text>\n"
                f"          </g>\n"
                f"          <g transform=\"translate(0,0),matrix(1 0 0 1 0 0)\">\n"
                f"            <g transform=\"translate(0,0),matrix(1 0 0 0.603015 0 0)\">\n"
                f"              <use href=\"url(#ideo_{path_id})\"/>\n"
                f"              <g clip-path=\"url(#ideo_{clip_path_id})\">\n"
                f"                <rect fill=\"url(#arm)\" height=\"{p_arm + 20}\" stroke=\"none\" width=\"12\" x=\"0\" y=\"0\"/>\n"
                #f"                <rect fill=\"url(#acen)\" height=\"{2}\" stroke=\"none\" width=\"12\" x=\"0\" y=\"{cen}\"/>\n"
                f"              </g>\n"
                f"            </g>\n"
                f"            <g transform=\"translate(0,0),matrix(1 0 0 0.603015 0 0)\"/>\n"
                f"          </g>\n"
                f"        </g>\n"
                )
                svg.write(chr_block)

    file_end = (
    "      </g>\n"
    "    </g>\n"
    "  </g>\n"  
    "</svg>\n"
    )

    svg.write(file_end)        
    
