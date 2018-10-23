'''
Author: Saideep Gona

This script is intended to populate the given database from local files
'''

from flask_app import db
from flask_app import ChIP_Meta, Peaks, Presets
import glob
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

# IO ************************************************************

pwd = os.getcwd()

if len(sys.argv) > 1:
    peak_dir = sys.argv[1]
    metadata_dir = sys.argv[2]
    metrics_dir = sys.argv[3]
else:
    peak_dir = pwd + "/temp_peaks/"
    metadata_path = pwd + "/tissue_types.pkl"
    metrics_dir =  pwd + "/static/images/"
# IO END *********************************************************

def find_duplicates(in_list):  
    unique = set(in_list)
    # print(unique)  
    for each in unique:  
        count = in_list.count(each)  
        if count > 1:  
            print ("duplicate: " + each + " " + count)

def histogram(field, data, metrics_dir):
    np_data = np.array(data)
    # print(np_data)
    # hist, bins = np.histogram(np_data, bins = 50)
    # print(hist,bins)
    plt.hist(data, color = 'blue', bins = 50)
    # # plt.plot(np.array([0,1,2,3]), color = 'blue')
    # plt.show()
    # print(x)
    savefile = metrics_dir + "/" + field + "_hist.png"
    print(savefile)
    plt.savefig(savefile)

# MAIN ***********************************************************

db.create_all()

# Read in metadata
metadata_dict = pickle.load(open(metadata_path, "rb"))
metadata_dict_ref = {}

find_duplicates(list(metadata_dict.keys()))
print(len(list(metadata_dict.keys())), " NUMBER OF STUDIES")
# sys.exit()

# print(metadata_files)
for m_f in list(metadata_dict.keys()):
    # print(m_f)

    tissue_obj = metadata_dict[m_f]["tissue"]
    if type(tissue_obj) == list:
        if len(tissue_obj) == 0:
            tissue = "NA"
        else:
            tissue = tissue_obj[-1]
    elif type(tissue_obj) == str:
        tissue = tissue_obj

    meta = ChIP_Meta(
        experiment_accession = m_f,
        tissue_types = tissue,
        transcription_factors = metadata_dict[m_f]["tf"]
    )

    metadata_dict_ref[m_f] = [tissue, metadata_dict[m_f]["tf"]]
    print(meta)
    db.session.add(meta)

db.session.commit()

# Read in all peak files

p_values = []
q_values = []
fold_enrichment = []

peak_files = glob.glob(peak_dir+"/*")
for p_f in peak_files:
    with open(p_f, "r") as pre_f:
        f = pre_f.readlines()[24:]
        for line in f:
            p_l = line.rstrip("\n").split("\t")
            if len(p_l) != 10:
                continue
            try:
                int(p_l[1])
            except:
                continue
            exp_acc = p_f.rstrip("_peaks.xls").split("/")[-1]
            if exp_acc not in metadata_dict_ref:
                print(exp_acc, " not in metadata!")
                continue
            peak = Peaks(
                experiment_accession = p_f.rstrip("_peaks.xls").split("/")[-1],
                tissue_types = metadata_dict_ref[exp_acc][0],
                transcription_factors = metadata_dict_ref[exp_acc][1],
                chrom = p_l[0],
                start = p_l[1],
                end = p_l[2],
                length = p_l[3],
                summit = p_l[4],
                pileup = p_l[5],
                log_p = p_l[6],
                fold_enrichment = p_l[7],
                log_q = p_l[8]
            )
            p_values.append(float(p_l[6]))
            q_values.append(float(p_l[8]))
            fold_enrichment.append(float(p_l[7]))
            print(peak)
            db.session.add(peak)
    db.session.commit()

print(len(p_values), max(p_values), min(p_values))
histogram("log_p", p_values, metrics_dir)
histogram("fold_enrichment", fold_enrichment, metrics_dir)


# END MAIN ********************************************************