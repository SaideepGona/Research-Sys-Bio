'''
Author: Saideep Gona

This script iterates through enhancer .fasta files in a given directory and creates a concatenated output file
which only contains the chrN:start-end information.

Enhancer data from the enhancer atlas should be converted from hg19 to hg38 using: https://genome.ucsc.edu/cgi-bin/hgLiftOver
'''

import os
import sys
from Bio import SeqIO

cwd = os.getcwd()

# IO ************************************************************

if len(sys.argv) > 1:

    raw_enhance_dir = sys.argv[1]
    out_file = sys.argv[2]

else:
    raw_enhance_dir = cwd + "/raw_enhancers/"
    out_file = "enhance_locations.txt"

# IO END ********************************************************


def pull_locs(filename):
    all_enhs = []
    # if filename.endswith("fasta"):
    #     f = SeqIO.parse(raw_enhance_dir + filename, "fasta")
    # else:
    #     f = open(raw_enhance_dir + filename, "r")
    # print(type(f))
    with open(raw_enhance_dir + filename, "r") as f:
        for line in f:
            if line[0] == '>':
                if line[1] == "H":
                    all_enhs.append(line.rstrip("\n").lstrip(">").split("|")[1].rstrip())
                else:
                    all_enhs.append(line.rstrip("\n").lstrip(">").split("_")[0].rstrip())
                keep_line = False
    return all_enhs

enh_raws = os.listdir(raw_enhance_dir)
full_enhs = []
for raw in enh_raws:
    full_enhs = full_enhs + pull_locs(raw)


if os.path.exists(out_file):
    os.remove(out_file)
with open(out_file, 'w') as out:
    for loc in full_enhs:
        out.write(loc + "\n")