'''
Author: Saideep Gona

Performs peak annotation for all .xls files in a folder (assumes they are for ChIP-Seq peaks or similar). 
Uses the "annotatePeaks.pl" script from HOMER.
'''

import glob
import os
import sys
import subprocess

# IO ************************************************************

if len(sys.argv) > 1:
    source_dir = sys.argv[1]
    out_dir = sys.argv[2]
else:
    source_dir = os.getcwd()
    out_dir = os.getcwd()

# IO END ********************************************************

all_xls_files = glob.glob(source_dir)

def annotate(exp):

    if exp[-3:] != ".xls":
        return
    annotate_out = out_dir + exp + ".anno"
    if os.path.exists(annotate_out):
        return
    make_file = "touch " + annotate_out
    os.system(make_file)
    annotation_file = open(annotate_out, 'w')
    annotate_process = subprocess.Popen(["annotatePeaks.pl", exp, "hg38"], stderr=None, stdout=annotation_file)
    output = annotate_process.communicate()
    annotation_file.close()


for xls in all_xls_files:
    annotate(xls)
