#!/usr/bin/env python3

from sys import argv

"""
Helper script for goldfleck pipeline

For a pair of simulated reads, this script logs the reference path,
percent of reads generated from said reference, and how many reads
that comes out to be.

Inputs: 3 positional arguments
- simrefs (str):    path to 'blend' containing the first 2 columns (ref path, percent)
- num_pairs (int):  sum total read pairs to generate
- outdir (str):     path to write resulting metadata file in

Output:
- outdir/metadata.csv

Author: Akash Pandit
Last Edited: Aug 29th, 2025
"""

simrefs   = str(argv[1])
num_pairs = int(argv[2])
outdir    = str(argv[3])

outlines = ['reference_path\tpercent_of_total_reads\tnumber_of_pairs_generated']

with open(simrefs, "r") as fp:
    lines = fp.read().strip().split('\n')
    delim = ',' if ',' in lines[0] else '\t'

    for line in lines:
        fp, percent = line.split(delim)
        fp, percent, numpairs = fp, int(percent), int(percent) * (num_pairs // 100)
        outlines.append(f"{fp}\t{percent}\t{numpairs}")

with open(f"{outdir}/blend-metadata.tsv", "w") as fp:
    fp.write('\n'.join(outlines) + '\n')
