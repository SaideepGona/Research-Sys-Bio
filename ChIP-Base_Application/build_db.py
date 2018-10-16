'''
Author: Saideep Gona

This script is intended to populate the given database from local files
'''

from flask_app import db
from flask_app import ChIP_Meta, Peaks, Presets
import glob
import sys
import pickle

# IO ************************************************************

if len(sys.argv) > 1:
    peak_dir = sys.argv[1]
    metadata_dir = sys.argv[2]
else:
    peak_dir = "/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/ChIP-Base_Application/peaks/"
    metadata_path = "/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/ChIP-Base_Application/tissue_types.pkl"

# IO END *********************************************************

def find_duplicates(in_list):  
    unique = set(in_list)
    # print(unique)  
    for each in unique:  
        count = in_list.count(each)  
        if count > 1:  
            print ("duplicate: " + each + " " + count)


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
            print(peak)
            db.session.add(peak)
    db.session.commit()


# END MAIN ********************************************************