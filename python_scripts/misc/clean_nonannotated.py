'''
Author: Saideep Gona

Deletes .tsv files of length 15 which have not been called for peaks
'''


import sys
import os
import glob

all_tsv_files = glob.glob("*.tsv")

for tsv_file in all_tsv_files:
    if len(tsv_file) == 15:
        if os.path.isfile(tsv_file[0:-4] + "_summits.bed") == False:
            os.remove(tsv_file)
            print(tsv_file)

