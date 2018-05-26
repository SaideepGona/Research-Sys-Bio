import pandas as pd
import glob
import subprocess
import os
import sys


if os.path.isfile('metadata_geo.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
    sys.stdout.flush()
    sys.exit()
else:
    full_meta_table = pd.read_table('metadata_geo.tsv')

if len(sys.argv) > 2:
    full_meta_table = pd.read_table(sys.argv[2])

completed_exps = [file_name[:-10] for file_name in glob.glob("*_peaks.xls")]   # Finds all completed experiments
completed_exps = completed_exps + [file_name[:-8] for file_name in glob.glob("*.skipped")]
print(completed_exps)
sys.stdout.flush()

completed_exps_bool = full_meta_table["series"].isin(completed_exps)  #Finds rows which are completed

incomplete_exps_bool = ~completed_exps_bool     # Finds rows which are not completed

incompleted_subset = full_meta_table[incomplete_exps_bool.values].reset_index(drop=True)    # Finds part of table which is not yet complete

incomplete_exp_count = incompleted_subset["series"].nunique() # Counts num incomplete experiments
partition_count = int(sys.argv[1])
partition_size = incomplete_exp_count // partition_count


current_partition_start = 0
experiments_included = 0
current_experiment = incompleted_subset.iloc[0]["series"]
current_partition = 0
for row in range(len(incompleted_subset)):
    if incompleted_subset.iloc[row]["series"] != current_experiment:
        experiments_included += 1
        current_experiment = incompleted_subset.iloc[row]["series"]
    if experiments_included == partition_size:
        incompleted_subset[current_partition_start:row].to_csv(str(current_partition+1) + "_inputdatageo.tsv", sep="\t", index = False)
        current_partition_start = row
        current_partition += 1
        experiments_included = 0
    if current_partition == partition_count - 1:
        incompleted_subset[current_partition_start:].to_csv(str(current_partition+1) + "_inputdatageo.tsv", sep="\t", index = False)
        break