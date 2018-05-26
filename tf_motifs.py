import pandas as pd
import numpy as np
import requests, json
import os
import math
import pickle

# COLLECTS TF MOTIFS(PWMs PREFERABLE) FOR ALL HUMAN TFS
# MASTER LIST COMES FROM http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19  <--NOT REALLY USED
#

cwd = os.getcwd()
jaspar_path = cwd + "/jaspar/"
hoco_path = cwd + "/hoco/pwm/"
all_tfs_path = cwd + "/tflist.txt"
pickle_filename = "tfs.pkl"
pickle_filename_all = "tfs_all.pkl"

skip_count = 0

def replace_zeros(matrix):
    new_matrix = np.zeros(matrix.shape)

    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            if matrix[row,col] == 0:
                new_matrix[row,col] = 0.01
            else:
                new_matrix[row,col] = matrix[row,col]

    return new_matrix

def parse_hoco(filename):
    '''
    Parse a single hocomoco file
        Inputs:
         filename: Path to target file
        Output:
         tf: name of transcription factor for associated weight matrix
         pwm: position weight matrix for transcription factor
    '''

    with open(filename) as f:
        content = []
        for line in f:
            content.append(line.rstrip("\n"))
    tf = content[0].lstrip(">").split("_")[0]
    split_content = []
    for line in content:
        split_content.append(line.split("\t"))
    pwm = np.array(split_content[1:])
    pwm = pwm.astype(float)
    print(pwm)
    return tf, pwm


def parse_all_hoco(pwm_dict, pwm_all_dict):
    '''
    Parse all hocomoco files in the hoco_path(global)
        Inputs:
         pwm_dict: dictionary of all tfs -> pwms
        Output:
         modified pwm_dict
    '''

    files = os.listdir(hoco_path)
    for fl in files:
        tf, pwm = parse_hoco(hoco_path + fl)

        if tf not in pwm_all_dict:
            pwm_all_dict[tf] = pwm

        if tf in pwm_dict:
            if type(pwm_dict[tf]) == bool:
                if pwm_dict[tf]:
                    pwm_dict[tf] = pwm


def parse_jaspar(filename):
    '''
    Parse a single jaspar file
        Inputs:
         filename: Path to target file
        Output:
         tf: name of transcription factor for associated weight matrix
         pwm: position weight matrix for transcription factor
    '''

    with open(filename) as f:
        content = []
        for line in f:
            content.append(line.rstrip("\n").split())
    # print(content)
    tf = content[0][1]
    pfm_t = np.array(content[1:])[:,2:-1]
    pfm = np.transpose(pfm_t)
    # print(tf)
    pwm = pfm_to_pwm(pfm)
    return tf, pwm

def pfm_to_pwm(pfm):
    '''
    Converts a position frequency matrix(pfm) to position weight matrix(pwm), assumes 0.25 background
        Inputs:
         pfm: Position frequency matrix to be converted (positions x bases)
        Output:
         pwm: position weight matrix for transcription factor post-conversion (positions x bases)
    '''

    totals = np.zeros(pfm.shape[0])
    pfm = pfm.astype(float)
    pfm = replace_zeros(pfm)

    for pos in range(pfm.shape[0]):
        # print(pfm[pos,:])
        totals[pos] = np.sum(pfm[pos,:])
    # print(totals)
    ppm = np.zeros(pfm.shape)

    for pos in range(pfm.shape[0]):
        for base in range(pfm.shape[1]):
            # print(totals[pos])
            # print(pfm[pos,base])
            ppm[pos,base] = pfm[pos,base]/totals[pos]

    pwm = np.zeros(pfm.shape)
    # print(ppm, "ppm")
    for pos in range(pfm.shape[0]):
        for base in range(pfm.shape[1]):
            pwm[pos,base] = math.log2(ppm[pos,base]/0.25)
    # print(pwm)
    return pwm, pwm

def parse_all_jaspar(pwm_dict, pwm_all_dict):
    '''
    Parse all jaspar files from jaspar path directory(global)
        Inputs:
         pwm_dict: dictionary of all tfs -> pwms
        Output:
         modified pwm_dict
         skip_count: number of skipped matrices due to numerical issues
    '''
    skip_count = 0
    files = os.listdir(jaspar_path)
    for fl in files:
        try:
            tf, pwm = parse_jaspar(jaspar_path + fl)
        except:
            skip_count += 1
            continue

        pwm_all_dict[tf] = pwm

        if tf in pwm_dict:
            if type(pwm_dict[tf]) == bool:
                if pwm_dict[tf]:
                    pwm_dict[tf] = pwm
    return skip_count


example_pfm = np.array([[1,2,3],[4,5,6],[1,2,3]])
print(pfm_to_pwm(example_pfm))

with open(all_tfs_path) as f:
    all_tfs = []
    for line in f:
        all_tfs.append(line.rstrip("\n"))

pwm_dict = {}
pwm_all_dict = {}
print(all_tfs)

for tf in all_tfs:
    pwm_dict[tf] = True

skips = parse_all_jaspar(pwm_dict,pwm_all_dict)
parse_all_hoco(pwm_dict,pwm_all_dict)

print(pwm_dict)
print(skips)

pwm_count = 0
for tf, val in pwm_dict.items():
    if type(val) != bool:
        pwm_count += 1

print(pwm_count, 'filled_entries')
print(len(pwm_all_dict))

# print(pwm_all_dict)


output = open(pickle_filename, 'wb')
pickle.dump(pwm_dict, output)
output.close()

output_all = open(pickle_filename_all, 'wb')
pickle.dump(pwm_all_dict, output_all)
output_all.close()

# print(pwm_dict)
