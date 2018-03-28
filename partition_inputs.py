import pandas as pd
import glob
import subprocess
import os
import sys

def download_link(url):
    """
    Download from link url. If there is an error along the way, tries again with continue invoked
    """

    try:
        print("downloading")
        subprocess.run(["wget", "-c", url])
    except:
        download_link(url)


files_txt_list = []
with open("files.txt") as f:                                # Open files.txt and pull out link to metadata.
    for line in f:
        files_txt_list.append(line.strip())

if os.path.isfile('metadata.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
    metadata = files_txt_list[0]
    download_link(metadata)
    full_meta_table = pd.read_table('metadata.tsv')
else:
    full_meta_table = pd.read_table('metadata.tsv')

completed_exps = [file_name[:-4] for file_name in glob.glob("*.tsv")]
print(completed_exps)
sys.stdout.flush()

completed_exps_bool = full_meta_table["Experiment accession"].isin(completed_exps)

incomplete_exps_bool = ~completed_exps_bool

incompleted_subset = full_meta_table[incomplete_exps_bool.values].reset_index(drop=True)

incomplete_exp_count = incompleted_subset["Experiment accession"].nunique()
partition_count = int(sys.argv[1])
partition_size = incomplete_exp_count // partition_count


current_partition_start = 0
experiments_included = 0
current_experiment = incompleted_subset.iloc[0]["Experiment accession"]
current_partition = 0
for row in range(len(incompleted_subset)):
    if incompleted_subset.iloc[row]["Experiment accession"] != current_experiment:
        experiments_included += 1
        current_experiment = incompleted_subset.iloc[row]["Experiment accession"]
    if experiments_included == partition_size:
        incompleted_subset[current_partition_start:row].to_csv(str(current_partition+1) + "_inputdata.tsv", sep="\t", index = False)
        current_partition_start = row
        current_partition += 1
        experiments_included = 0
    if current_partition == partition_count - 1:
        incompleted_subset[current_partition_start:].to_csv(str(current_partition+1) + "_inputdata.tsv", sep="\t", index = False)
        break




