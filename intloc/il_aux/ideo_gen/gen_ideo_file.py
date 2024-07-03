# Copyright 2022 Simon Magin. 
# Licensed under the BSD-2-clause License (https://opensource.org/licenses/BSD-2-Clause).
# This file may not be copied, modified, or distributed except according to those terms.

import sys
import os

svg = sys.argv[1]

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
    "  </defs>"
    "    <g transform=\"translate(0,0),matrix(1 0 0 1 0 0)\">\n"
    "    <g transform=\"translate(0,0),matrix(1 0 0 1 0 0)\">\n"
    "      <g transform=\"translate(20,20),matrix(1 0 0 1 0 16)\">\n"
    )

    svg.write(preamble)

    n = 22

    for i in range(n):
        transl = i * 24
        i = i+1
        clip_path_id = i*2
        path_id = clip_path_id -1
        chrom = i
        chr_block = (
        f"    <clipPath id=\"ideo_{clip_path_id}\">"
        f"      <path id=\"ideo_{path_id}\" d=\"M0,5.05 C0,2.5 3,0 6,0 C9,0 12,2.5 12,5 L12,160 C12,162 9,165 6 165 C3,165 0,162 0,160 Z M0,170 C0,167 3,165 6,165 C9,165 12,167 12,170 L12,332 C12,335 9,337 6 337 C3,337 0,335 0,332 Z\" fill=\"rgb(240,240,255)\" stroke=\"rgb(0,0,0)\" stroke-width=\"1.000px\"/>\n"
        f"    </clipPath>\n"  
        f"        <g transform=\"translate({transl},0),matrix(1 0 0 1 0 0)\">\n"
        f"          <g fill=\"rgb(32,32,32)\" fill-opacity=\"1\" transform=\"translate(3,-16),matrix(1 0 0 1 0 0)\">\n"
        f"            <text font-family=\"Courier\" font-size=\"12.000000px\" transform=\"matrix(1 0 0 1 0 8)\" x=\"0\" y=\"0\">{chrom}</text>\n"
        f"          </g>\n"
        f"          <g transform=\"translate(0,0),matrix(1 0 0 1 0 0)\">\n"
        f"           <g transform=\"translate(0,0),matrix(1 0 0 0.603015 0 0)\">\n"
        f"              <use href=\"url(#ideo_{path_id})\"/>\n"
        f"              <g clip-path=\"url(#ideo_{clip_path_id})\">\n"
        f"                <rect fill=\"url(#arm)\" height=\"357\" stroke=\"none\" width=\"12\" x=\"0\" y=\"0\"/>\n"
        f"                <rect fill=\"url(#acen)\" height=\"20\" stroke=\"none\" width=\"12\" x=\"0\" y=\"155\"/>\n"
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
    
