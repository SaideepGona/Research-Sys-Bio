'''
Author: Saideep Gona

This script is intended to create quickly searchable bed files from chip-seq data, as well as organize some associated 
dnase and enhancer files.
'''

# This script is intended to create quickly searchable bed files from chip-seq data

import os
import sys
import pandas as pd
import numpy as np
import glob

pwd = os.getcwd()

# IO ************************************************************

if len(sys.argv) > 1:
    bed_dir = sys.argv[1]
    encode_data_dir = pd.read_table(sys.argv[2])
    dnase_dir = sys.argv[3]
    enhancer_dir = sys.argv[4]
else:
    bed_dir = pwd + "/bed_files"
    encode_data = pd.read_table(pwd + '/metadata.tsv')
    dnase_dir = bed_dir +"/AA_dnase/"
    enhancer_dir = bed_dir +"/AA_enhancers/"

# IO END *********************************************************

def copy_bed(source, destination):

    command = "\cp "+source+" "+destination
    os.system(command)

def remove(path):
    if not os.path.exists(path):
        return
    else:
        os.remove(path)

def annotation_to_bed(anno, distance_threshold_upstream, distance_threshold_downstream):

    annotation_df = pd.read_csv(anno, sep='\t')

    # Insert filtering criteria below
    distance_df = (annotation_df[(annotation_df["Distance to TSS"] > distance_threshold_upstream) \
                        & (annotation_df["Distance to TSS"] < distance_threshold_downstream) \
                        ])

    all_genes = []
    for index, row in distance_df.iterrows():
        all_genes.append(str(row["Gene Name"]))

    write_list = [self.transcription_factor] + all_genes

    out_file = open(self.accession_id + ".tfgenes", 'w')
    out_file.write("\n".join(write_list) + "\n")
    out_file.close()

with open(pwd + '/tflist.txt') as tfs_file:                 # Make a list of all possible tf names
    tfs_pre = tfs_file.read()
    tfs = tfs_pre.split("\n")

tf_counts = {}                                              # Initialize tf counts dictionary
for tf in tfs:
    tf_counts[tf] = 0

os.makedirs(bed_dir, exist_ok=True)                         # Make directory for each
for tf in tfs:
    os.makedirs(bed_dir+"/"+tf, exist_ok=True)

for row in encode_data.itertuples():                        # Copy bed files into the correct corresponding directories
    # print(row[17])
    tf_raw = row[17]
    accession = row[4]
    tf = tf_raw.rstrip('-human')

    if tf in tf_counts:
        source = bed_dir + "/" + accession + "_peaks.xls"
        dest = bed_dir + "/" + tf + "/"
        target = dest + accession + "_peaks.xls"
        if os.path.exists(target):                          # If already exists in place dont actually do anything
            os.remove(target)
            continue
        copy_bed(source, dest)
        tf_counts[tf] += 1

# print(tfs)

for tf in tfs:                                               # Merge .bed files in each tf folder and sort the output
    if tf == None:
        continue
    if tf == " ":
        continue
    if tf == "":
        continue

    print(tf)

    remove(bed_dir+"/"+tf+"/"+tf+"_sorted.bed")
    remove(bed_dir+"/"+tf+"/"+tf+"_sorted.bed.final")

    merge_beds_full = glob.glob(bed_dir+"/"+tf+"/*")
    if len(merge_beds_full) == 0:
        continue
    # print(tf, "tf")
    # merge_beds_full = [(bed_dir+"/"+x) for x in merge_beds]
    # print(merge_beds_full)

    sorted_file = bed_dir + "/" + tf + "/" + tf + "_sorted.bed"
    input_files = " ".join(merge_beds_full)
    # print(input_files)
    cat_command = "cat " + input_files + " > " + sorted_file
    os.system(cat_command)

    sort_command = "bedtools sort -i " + sorted_file + " > " + sorted_file+".final"
    os.system(sort_command)

# Keep enhancer regions ready
print("enhancer")
sorted_file = enhancer_dir + "enhancer_sorted.bed"
remove(sorted_file)
remove(sorted_file+".final")

sort_command = "bedtools sort -i " + "full_enhancer_locations_hg38.bed" + " > " + sorted_file+".final"
os.system(sort_command)


# Keep dnase footprints ready
print("dnase")
sorted_file = dnase_dir + "dnase_sorted.bed"
remove(sorted_file)
remove(sorted_file+".final")

dnase_files = glob.glob(dnase_dir+"*")
input_dnase = " ".join(dnase_files)

cat_command = "cat " + input_dnase + " > " + sorted_file
os.system(cat_command)

sort_command = "bedtools sort -i " + sorted_file + " > " + sorted_file+".final"
os.system(sort_command)

annotate_command = "annotatePeaks.pl "+sorted_file+".final " + "hg38 > "+dnase_dir+"dnase_annotations.txt"
os.system(annotate_command)

unique = 0
for tf, count in tf_counts.items():
    if count > 0:
        unique += 1
print(unique)










