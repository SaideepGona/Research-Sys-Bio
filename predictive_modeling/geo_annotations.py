import os
import pandas as pd
import numpy as np
import subprocess
import urllib
import sys
import requests, json
import glob
import time

# Dependencies:
# python libraries above
# wget
# bowtie2
# samtools
# macs2

# Start with files.txt file

# TO-DO
# TODO Detect whether given experiment already has .bam or peaks
# TODO Create alternative pipelines as needed
# TODO Use Hg38 reference genome

distance_threshold_upstream = -1000                         #Permissable peak distances from transcription start site for inclusion
distance_threshold_downstream = 100
time = str(time.time())
scratch_path_part = "/scratch/sgona/"
scratch_path = scratch_path_part + time +"/"
# scratch_path = ""
subprocess.run(["mkdir", "-p", scratch_path_part])
subprocess.run(["mkdir", "-p", scratch_path+"/tmp/"])
# scratch_path = "/home/saideep/Documents/GitHub_Repos/Saideep/MSCB_Sem1/Research/Research-Sys-Bio/"+time+"/"
# subprocess.run(["mkdir", scratch_path])


class Experiment:

    def __init__(self, accession_id, transcription_factor, noncontrols, controls):

        # For GEO, the "series" accesssion is the equivalent of the experiment-wide accession
        # Transcription factor should be assessed in read-in file
        # I think GEO processing will start directly from fastq b/c reference genome is rarely mentioned

        self.accession_id = accession_id
        self.transcription_factor = transcription_factor
        self.noncontrols = noncontrols
        self.controls = controls
        self.file_state = "fastq"
        self.processed_file = None

    # def assess_experiment(self):
    #     #Checks if .bam/peaks are present and stores this as a property)
    #     exp_json = get_request(self.accession_id)                           # Gets json for given experiment and converts to dict
    #     files = exp_json['files']                                           # Pulls out file list from expreiment json
    #     files_acc = [(file_string[7:])[:-1] for file_string in files]
    #     file_props = [find_file_props(acc) for acc in files_acc]
    #     files = self.determine_file_state(file_props)

    # def determine_file_state(self, file_props):
    #     """
    #     Determines the current level of processing for the experiment(fastq,bam,bed)
    #     Also returns a list of accessions for all files of the most-processed file type
    #     """
    #     def is_fully_processed(self, file_props):
    #         #first determine biological replicates
    #         largest_group = []

    #         for ind_file in file_props:
    #             if ind_file[3] == None:
    #                 continue
    #             if len(ind_file[3]) > len(largest_group):
    #                 largest_group = ind_file[3]

    #         target_properties = ["bed narrowPeak", "peaks", largest_group]

    #         for ind_file in file_props:
    #             if ind_file[1:] == target_properties:
    #                 self.processed_file = ind_file[0]
    #                 return True

    #     fastq_files = []
    #     bam_files = []
    #     peak_files = []

    #     if is_fully_processed(self, file_props):        #Handles the fully processed case
    #         self.file_state = "bed"


    #     # [accession, file_type, output_type, br] <- file preps is a list of these lists

    #     for single_file in file_props:
    #         file_type = single_file[1]
    #         if single_file[1] == "fastq" and len(single_file[3]) == 1:
    #             fastq_files.append(single_file[0])
    #         elif single_file[1] == "bam" and single_file[2] == "alignments":
    #             replicate = Replicate(single_file[0], False)
    #             bam_files.append(replicate)

    #     if len(bam_files) > 0:
    #         self.file_state = "bams"
    #         self.bam_noncontrol_files = bam_files
    #         self.bam_control_files = self.find_bam_controls()
    #         return bam_files

    #     elif len(fastq_files) > 0:
    #         self.file_state = "fastq"
    #         return fastq_files

    # def find_bam_controls(self):

    #     print("find_bam_controls")
    #     sys.stdout.flush()

    #     # Finds all control .bam files for the given experiment

    #     response = get_request(self.accession_id)
    #     controls = response["possible_controls"]
    #     controls_acc = [(control_acc[13:])[:-1] for control_acc in controls]
    #     all_control_bams = []
    #     all_control_replicates = []
    #     print(controls_acc, "control accessions, can have files within them")
    #     for control_exp in [controls_acc[0]]:
    #         control_bams = self.find_bams_in_controlexps(control_exp)
    #         print(control_bams
    #         all_control_bams = all_control_bams + control_bams
    #         all_control_bams = list(set(all_control_bams))

    #     for single_rep in all_control_bams:
    #         all_control_replicates.append(Replicate(single_rep, True))

    #     return all_control_replicates


    # def find_bams_in_controlexps(self, control_acc):

    #     # Given a control experiment accession -> outputs all bam control files
    #     exp_json = get_request(control_acc)                           # Gets json for given experiment and converts to dict
    #     files = exp_json['files']                                           # Pulls out file list from experiment json
    #     files_acc = [(file_string[7:])[:-1] for file_string in files]
    #     file_props = [find_file_props(acc) for acc in files_acc]

    #     target_properties = ["bam", "alignments"]

    #     control_bams = []

    #     for ind_file in file_props:
    #         if ind_file[1:3] == target_properties:
    #             control_bams.append(ind_file[0])
    #     print(control_bams)
    #     return control_bams

    def process_experiment(self):

        '''
        Given the "state" of the experiment, performs steps until .bed output for the experiment
        '''
        self.file_state = "fastq"
        # sys.stdout.flush()
        # # Handles if file has already been processed
        # if self.file_state == "bed":
        #     download_link(accession_to_url(self.processed_file,".bed.gz"))
        #     os.rename(scratch_path + self.processed_file + ".bed.gz", scratch_path + self.accession_id + ".bed.gz")
        #     subprocess.run(["gunzip", scratch_path + self.accession_id + ".bed.gz"])
        #     #sys.stdout.flush()
        #     return

        # elif self.file_state == "bams":

        #     merged_noncontrols = self.create_merged_bam_no_alignment(False)
        #     merged_controls = self.create_merged_bam_no_alignment(True)
        #     print(merged_controls, merged_noncontrols, "BOTH MERGED BAM FILES")
        #     sys.stdout.flush()
        #     self.create_macs_peak(merged_noncontrols, merged_controls)
        #     #sys.stdout.flush()

        if self.file_state == "fastq":
            merged_noncontrols = self.create_merged_bam(False)
            merged_controls = self.create_merged_bam(True)
            if merged_noncontrols == None or merged_noncontrols == None:
                return
            self.create_macs_peak(merged_noncontrols, merged_controls)

    # def create_merged_bam_no_alignment(self, is_control):

    # # Given if we are working with reps or controls, takes the list of files and processes from fastq to merged bam
    # # depending on the current level of processing
    #     print("MERGING BAMS-no alignment")
    #     sys.stdout.flush()
    #     if is_control:
    #         current_reps = self.bam_control_files
    #         merged_filename = "_merged_controls.bam"
    #     else:
    #         current_reps = self.bam_noncontrol_files
    #         merged_filename = "_merged_noncontrols.bam"

    #     rep_bam_files = []
    #     for rep in current_reps:
    #         bam_file = rep.download_only_bam()
    #         rep_bam_files.append(bam_file)              # List of (accession + .bam) filenames to be merged

    #     merged_file = scratch_path + self.accession_id + merged_filename

    #     if len(rep_bam_files) == 1:
    #         merged_file = rep_bam_files[0]
    #     elif len(rep_bam_files) > 1:
    #         subprocess.run(["stdbuf", "-o", "500MB","samtools", "merge", merged_file]+rep_bam_files)

    #     for rep_file in rep_bam_files:
    #         delete_file(rep_file)
    #         delete_file(rep_file[:-4]+".sam")
    #         delete_file(rep_file[:-4]+".fastq")
    #         delete_file(rep_file[:-4]+".place")

    #     return merged_file


    def create_merged_bam(self, is_control):

        # Given if we are working with reps or controls, takes the list of files and processes from fastq to merged bam
        # depending on the current level of processing
        print('merging bams')
        if is_control:
            current_reps = self.controls
            merged_filename = "_merged_controls.bam"
            print(self.controls)
        else:
            current_reps = self.noncontrols
            merged_filename = "_merged_noncontrols.bam"
            print(self.noncontrols)

        sys.stdout.flush()

        rep_bam_files = []
        for rep in current_reps:
            bam_file = rep.download_to_bam()
            if bam_file == None:
                subprocess.run(["touch", scratch_path + self.accession_id + ".tsv"])
                return
            rep_bam_files.append(bam_file)              # List of (accession + .bam) filenames to be merged

        merged_file = scratch_path + self.accession_id + merged_filename

        if len(rep_bam_files) == 1:
            merged_file = rep_bam_files[0]
        elif len(rep_bam_files) > 1:
            subprocess.run(["samtools", "merge", merged_file]+rep_bam_files)

        for rep_file in rep_bam_files:
            if len(rep_bam_files) > 1:
                delete_file(rep_file[:-4]+".bam")
            delete_file(rep_file[:-4]+".sam")
            delete_file(rep_file[:-4]+".fastq")
            delete_file(rep_file[:-4]+".place")

        return merged_file

    def create_macs_peak(self, merged_noncontrols, merged_controls):
        #print(glob.glob(scratch_path + "*"))
        print("DIRECTLY PRIOR TO PEAK CALLING FILE STATE IS", self.file_state)
        print(merged_controls,merged_noncontrols)
        sys.stdout.flush()

        subprocess.run(["macs2", "callpeak", "--tempdir", scratch_path + "/tmp/", "-t", merged_noncontrols, "-c", merged_controls, "-n", self.accession_id, "--outdir", scratch_path]) # Call Macs Peak of the merged files
        subprocess.run(["cp", scratch_path + self.accession_id + "_summits.bed", "/home/sgona/encode_scraper"])

        delete_file(merged_noncontrols)
        delete_file(merged_controls)

    def annotate_experiment(self):

       # print(glob.glob(scratch_path + "*"))
        print("DIRECTLY PRIOR TO ANNOTATION FILE STATE IS", self.file_state)
        sys.stdout.flush()

        if self.file_state == "bed":
            experiment_bed_file = scratch_path + self.accession_id + ".bed"
        else:
            experiment_bed_file = scratch_path + self.accession_id + "_summits.bed"

        annotate_out = self.accession_id + ".tsv"
        annotation_file = open(annotate_out, 'w')
        annotate_process = subprocess.Popen(["annotatePeaks.pl", experiment_bed_file, "hg38"], stderr=None, stdout=annotation_file)
        output = annotate_process.communicate()
        annotation_file.close()

    def process_annotation(self):

        #print(glob.glob(scratch_path + "*"))
        print("DIRECTLY PRIOR TO ANNOTATION PROCESSING FILE STATE IS", self.file_state)
        sys.stdout.flush()

        annotation_df = pd.read_csv(self.accession_id + ".tsv", sep='\t')
        threshold_file = open(scratch_path + str(distance_threshold_upstream)+"_"+str(distance_threshold_downstream)+".tableprobs", 'w') #Make tableprops file
        threshold_file.close()
        # Insert filtering criteria below
        distance_df = (annotation_df[(annotation_df["Distance to TSS"] > distance_threshold_upstream) \
                            & (annotation_df["Distance to TSS"] < distance_threshold_downstream) \
                            ])

        all_genes = []
        for index, row in distance_df.iterrows():
            all_genes.append(str(row["Gene Name"]))

        write_list = [self.transcription_factor] + all_genes

        out_file = open(self.accession_id + ".tfgenes", 'w')
        sys.stdout.flush()
        out_file.write("\n".join(write_list) + "\n")
        out_file.close()

class Replicate:

    def __init__(self, accession_id, is_control, sra_path):

        #Replicate accession IDs are GEO sample accessions
        self.accession_id = accession_id
        self.is_control = is_control
        self.sra_path = sra_path

    def download_to_bam(self):
        """
        Downloads a fastq for a given accession ID, and processes it into a corresponding .bam file
        """

        subprocess.run(["touch", scratch_path + self.accession_id + ".place"])              # Creates a .place file as a placeholder
        download_link(self.sra_path)
        if not os.path.exists(scratch_path + self.sra_path+".fastq"):
            return
        os.rename(scratch_path + self.sra_path + ".fastq", scratch_path + self.accession_id + ".fastq")
        download_file = scratch_path + self.accession_id + ".fastq"
        # Align .fastq and output .sam; .sam is piped into samtools view to output .bam

        bam_out = scratch_path + self.accession_id + ".bam"
        bam_file = open(bam_out, 'w')
        alignment_process = subprocess.Popen(["bowtie2", "-x ", reference_genome, "-U", download_file], bufsize = 500*10**6, stderr=None, stdout=subprocess.PIPE)
        sam_to_bam_process = subprocess.Popen(["stdbuf", "-o", "1000MB", "samtools", "view", "-bSu", "-"], stderr=None, stdin=alignment_process.stdout, stdout=bam_file)
        alignment_process.stdout.close()
        output = sam_to_bam_process.communicate()
        alignment_process.wait()

        delete_file(download_file)                                          # Convert .sam to .bam

        quick_check(bam_out)
        print(self.accession_id, 'accession id')
        print(bam_out)
        sys.stdout.flush()
        return bam_out

    def download_only_bam(self):

        download_link(accession_to_url(self.accession_id, ".bam"))
        return scratch_path + self.accession_id + ".bam"

def get_request(accession):

    headers = {'accept': 'application/json'}
    url = "http://www.encodeproject.org/biosample/" + accession + "/?frame=object"
    response = requests.get(url, headers = headers)
    response_dict = response.json()

    return response_dict

def find_file_props(accession):

    response = get_request(accession)
    file_type = response.get("file_type")
    br = response.get("biological_replicates")
    output_type = response.get("output_type")

    return [accession, file_type, output_type, br]

def quick_check(file):
    """
    Runs the samtools quickckeck command on a given input file to check for truncation. If there is a problem,
    will cancel the job and write the filename to stdout
    """
    check_out = subprocess.run(["samtools", "quickcheck", "-v", file])

    if check_out.returncode != 0:
        print("quickcheck failed")
        sys.stderr.write("Truncated File:" + file)
        sys.stdout.flush()
        sys.exit()

def download_link(url):
    """
    Download from link url. If there is an error along the way, tries again with continue invoked
    """


    #subprocess.run(["wget", "-nv","-c", url, "-P", scratch_path])
    print('downloading')
    subprocess.run(["fastq-dump", "--outdir", scratch_path, url])             #use fastq-dump for downloading from GEO


def delete_file(filename):
    """
    # Delete a file if it exists. Otherwise do nothing
    """

    if os.path.isfile(filename) == True:
        os.remove(filename)
    else:
        return

# def accession_to_url(accession, file_tag):
#     """
#     Converts a file accession to a relevant download URL
#     """

#     url = 'https://www.encodeproject.org/files/' + accession + '/@@download/' + accession + file_tag
#     return url


def experiment_exists(experiment_ID, experiment):
    """
    Checks if an experiment has already been processed or is being processed
    """

    def extension_exists(root):
        rep_bool = (os.path.isfile(root + ".fastq.gz")
            or os.path.isfile(root + ".fastq")
            or os.path.isfile(root + ".sam")
            or os.path.isfile(root + ".bam")
            or os.path.isfile(root + ".place")
        )
        if rep_bool:
            print("exists")
            return True
        return False


    for noncontrol in experiment.noncontrols:
        if extension_exists(scratch_path + noncontrol.accession_id):
            return True
    for control in experiment.controls:
        if extension_exists(scratch_path + control.accession_id):
            return True

    exp_bool = (os.path.isfile(scratch_path+experiment_ID + "-peaks.tsv")
        or os.path.isfile(scratch_path+experiment_ID + ".place")
        or os.path.isfile(scratch_path+experiment_ID + ".bed.gz")
        or os.path.isfile(scratch_path+experiment_ID + ".bed")
        or os.path.isfile(scratch_path+experiment_ID + "-narrowpeak.bed")
        or os.path.isfile(scratch_path+experiment_ID + "-summits.bed")
        or os.path.isfile(scratch_path+experiment_ID + "-log.txt")
        or os.path.isfile(scratch_path+experiment_ID + "-model.pdf")
        or os.path.isfile(scratch_path+experiment_ID + "_peaks.narrowpeak")
        or os.path.isfile(scratch_path+experiment_ID + "_peaks.xls")
        or os.path.isfile(scratch_path+experiment_ID + "_summits.bed")
        or os.path.isfile(scratch_path+experiment_ID + "_model.r")
    )
    if exp_bool:
        print(experiment_ID, "exists")
        #sys.stdout.flush()
        return True
    return False

def write_bad_exp(exp_accession,note):

    """
    Logs bad experiments with some odd properties to be revisited later
    """

    f = open("bad_experiments.txt", "a")
    f.write(str(exp_accession) + " - " + note)
    f.close()


# Main Code: Requires a files.txt file which can be downloaded from the Encode website after selecting by field in the data Matrix


# DOWNLOAD METADATA AND READ IN AS A DATAFRAME **********************************************************************

input_file = sys.argv[1]
print(input_file)
meta_table = pd.read_table(input_file)
print(meta_table)
sys.stdout.flush()

#files_txt_list = []
#with open("files.txt") as f:                                # Open files.txt and pull out link to metadata.
#    for line in f:
#        files_txt_list.append(line.strip())

#if os.path.isfile('metadata.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
#    metadata = files_txt_list[0]
#    download_link(metadata)
#    meta_table = pd.read_table('metadata.tsv')
#else:
#    meta_table = pd.read_table('metadata.tsv')               # If metadata.tsv exists, read as dataframe


# ******************************************************************************************************************
# DOWNLOAD HG38 REFERENCE GENOME ***********************************************************************************

global reference_genome
reference_genome = "hg38ind"
ref_download_filename = "hg38.chromFa.tar.gz"
if os.path.isfile(ref_download_filename) == False:
    subprocess.run(["touch", reference_genome + ".place"])
    ref_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/" + ref_download_filename
    download_link(ref_url)
    subprocess.run(["gunzip", ref_download_filename])
    subprocess.run(["tar", "-xkf", ref_download_filename[0:-3]])
    subprocess.run(["cat", "chroms/*.fa", ">", "reference.fa"])
    subprocess.run(["bowtie2-build", "reference.fa", reference_genome])         #Builds reference genome

# *****************************************************************************************************************
# GENERATE EXPERIMENT TO TRIAL MATCHING **************************************************************************

value_counts = meta_table['series'].value_counts() # Finds all unique experiments and number of associated files for each
#print(type(valueCounts), len(valueCounts), "Experiments Conducted")

# Generate mapping of all accession IDs per experiment to each transcription factor
global accession_length
accession_length = 11            # Length of file accession IDs (can change over time?)

experiments_dict = {}             # For each key(Experiment Accession ID), there is a corresponding Experiment Object

for exp, val in value_counts.iteritems():

    experiments_dict[exp] = Experiment(exp, None, [], [])                                  # Add experiment object

    for key, replicate in meta_table.loc[meta_table['series'] == exp, 'iid'].iteritems():  # Add replicates and controls
        # print(key, replicate)
        sys.stdout.flush()
        sra_paths_str = meta_table['sra_path'][key]
        if type(sra_paths_str) != str:
            continue
        print(sra_paths_str)
        sys.stdout.flush()
        sra_paths = sra_paths_str.split(',')
        print(meta_table['is_control'][key], "whyyy")
        sys.stdout.flush()                                        # to experiment object
        if meta_table['is_control'][key] == False:
            print('true')
            for sra_path in sra_paths:
                experiments_dict[exp].noncontrols.append(Replicate(replicate, False, sra_path))
        elif meta_table['is_control'][key] == True:
            for sra_path in sra_paths:
                experiments_dict[exp].controls.append(Replicate(replicate, True, sra_path))

    # for key, control in meta_table.loc[meta_table['Experiment accession'] == exp, 'Controlled by'].iteritems():

    #     if type(control) is not str:
    #         continue

    #     if len(control) > 2*accession_length:                               # Some controls are strangely organized
    #         all_controls = control.split(',')
    #         for indiv_control in all_controls:
    #             accession_only = indiv_control[-(accession_length+1):-1]
    #             experiments_dict[exp].controls.append(Replicate(accession_only, is_control = True))
    #     else:
    #         accession_only = control[-(accession_length+1):-1]
    #         control_download = accession_to_url(accession_only, "fastq.gz")
    #         experiments_dict[exp].controls.append(Replicate(accession_only, is_control = True))

# print(experiments_dict)
# print("number of studies", len(experiments_dict))
# print(meta_table["Experiment target"].value_counts())
# *****************************************************************************************************************
# CREATE MACS PEAK AND ANNOTATIONS FOR EACH EXPERIMENT ****************************************************************************

bad_experiments = []
for exp, experiment in experiments_dict.items():
    print(exp)
    print(experiment.controls)
    print(experiment.noncontrols)
    if len(experiment.controls) == 0:
        print('no controls')                                   # Check if replicates or controls are empty
        write_bad_exp(exp, 'no control')
        continue
    if len(experiment.noncontrols) == 0:
        print('no reps')
        write_bad_exp(exp, 'no replicates')
        continue
    sys.stdout.flush()
    subprocess.run(["touch", scratch_path+experiment.accession_id + ".place"])
    print(experiment.accession_id)
    sys.stdout.flush()
    # experiment.assess_experiment()
    # process_exp = experiment.process_experiment()
    # if process_exp == None:
        # continue
    experiment.annotate_experiment()
    #experiment.process_annotation()

