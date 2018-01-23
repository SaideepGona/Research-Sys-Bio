import wget
import os
import pandas as pd
import numpy as np
import subprocess
import urllib
import sys
import requests, json

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

class Pipeline:



class Experiment:

    def __init__(self, accession_id, transcription_factor, noncontrols, controls):

        self.accession_id = accession_id
        self.transcription_factor = transcription_factor
        self.noncontrols = noncontrols
        self.controls = controls

    def assess_experiment(self):
        #Checks if .bam/peaks are present and stores this as a property
        exp_json = get_request(self.accession_id)
        files = exp_json['files']
        files_acc = [(file_string[7:])[:-1] for file_string in files]
        file_types = [find_file_type(file_acc) for file_acc in files]

        files = determine_file_state(file_types)


    def determine_file_state(self, files, file_types):
        """
        Determines the current level of processing for the experiment(fastq,bam,peak)
        Also returns a list of accessions for all files of the most-processed file type
        """

        fastq_files = []
        bam_files = []
        peak_files = []

        for x in range(len(file_types)):
            if file_types[x] == "fastq":
                fastq_files.append(files[x])
            elif file_types[x] == "bam":
                fastq_files.append(files[x])
            elif file_types[x] == "bed":
                peak_files.append(files[x])

        if len(peak_files) > 0:
            self.file_state = "peaks"
            return peak_files
        elif len(bam_files) > 0:
            self.file_state = "bams"
            return bam_files
        elif len(fastq_files) > 0:
            self.file_state = "fastqs"
            return fastq_files

    def process_experiment(self):

        merged_noncontrols = self.create_merged_bam(self.accession_id, False)
        merged_controls = self.create_merged_bam(self.accession_id, True)
        self.create_macs_peak(merged_noncontrols, merged_controls)

    def create_merged_bam(self, experiment_ID, is_control):

        if is_control:
            current_reps = self.controls
            merged_filename = "_merged_controls.bam"
        else:
            current_reps = self.noncontrols
            merged_filename = "_merged_noncontrols.bam"

        rep_bam_files = []
        for rep in current_reps:
            bam_file = rep.download_to_bam()
            rep_bam_files.append(bam_file)

        merged_file = experiment_ID + merged_filename

        if len(rep_bam_files) == 1:
            merged_file = rep_bam_files[0]
        elif len(rep_bam_files) > 1:
            subprocess.run(["samtools", "merge", merged_file]+rep_bam_files)

        for rep_file in rep_bam_files:
            delete_file(rep_file)
            delete_file(rep_file[:-4]+".sam")
            delete_file(rep_file[:-4]+".fastq")
            delete_file(rep_file[:-4]+".place")

        return merged_file

    def create_macs_peak(self, merged_noncontrols, merged_controls):

        try:
            subprocess.run(["macs2", "callpeak", "-t", merged_noncontrols, "-c", merged_controls, "-n", self.accession_id]) # Call Macs Peak of the merged files
        except:
            sys.stdout.write("Problem running MACS2")
            sys.exit()

        delete_file(merged_noncontrols)
        delete_file(merged_controls)


class Replicate:

    def __init__(self, accession_id, is_control):

        self.accession_id = accession_id
        self.is_control = is_control

    def download_to_bam(self):
        """
        Downloads a fastq for a given accession ID, and processes it into a corresponding .bam file
        """

        subprocess.run(["touch", self.accession_id + ".place"])              # Creates a .place file as a placeholder
        download_link(accesion_to_url(self.accession_id))
        download_file = self.accession_id + ".fastq.gz"
        print("downloaded file")
        subprocess.run(["gunzip", download_file])                       # Unzip the .fastq.gz to a fastq
        print("unzipped")
        delete_file(download_file)
        # Align .fastq and output .sam; .sam is piped into samtools view to output .bam

        bam_out = self.accession_id + ".bam"
        bam_file = open(bam_out, 'w')
        alignment_process = subprocess.Popen(["bowtie2", "-x ", reference_genome, "-U", download_file[0:-3]], bufsize = 500*10**6, stderr=None, stdout=subprocess.PIPE)
        sam_to_bam_process = subprocess.Popen(["samtools", "view", "-bSu", "-"], stderr=None, stdin=alignment_process.stdout, stdout=bam_file)
        alignment_process.stdout.close()
        output = sam_to_bam_process.communicate()
        alignment_process.wait()

        delete_file(download_file[0:-3])                                          # Convert .sam to .bam

        quick_check(bam_out)

        return bam_out

def get_request(accession):

    headers = {'accept': 'application/json'}
    url = "http://www.encodeproject.org/biosample/" + accession + "/?frame=object"
    response = requests.get(url, headers = headers)
    response_dict = response.json()

    return response_dict

def find_file_type(accession):

    response = get_request(accession)
    file_type = response["file_format"]

    return file_type

def quick_check(file):
    """
    Runs the samtools quickckeck command on a given input file to check for truncation. If there is a problem,
    will cancel the job and write the filename to stdout
    """
    check_out = subprocess.run(["samtools", "quickcheck", "-v", file])
    print(check_out)

    if check_out.returncode != 0:
        print("quickcheck failed")
        sys.stderr.write("Truncated File:" + file)
        sys.exit()

def download_link(url):
    """
    Download from link url. If there is an error along the way, tries again with continue invoked
    """

    try:
        print("downloading")
        subprocess.run(["wget", "-c", url])
    except:
        download_link(url)

def delete_file(filename):
    """
    # Delete a file if it exists. Otherwise do nothing
    """

    if os.path.isfile(filename) == True:
        os.remove(filename)
        sys.stdout.write(filename + "deleted")
    else:
        return

def accesion_to_url(accession):
    """
    Converts a file accession to a relevant download URL
    """

    url = 'https://www.encodeproject.org/files/' + accession + '/@@download/' + accession + '.fastq.gz'
    return url


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
        print(rep_bool)
        return False


    for noncontrol in experiment.noncontrols:
        if extension_exists(noncontrol.accession_id):
            return True
    for control in experiment.controls:
        if extension_exists(control.accession_id):
            return True

    exp_bool = (os.path.isfile(experiment_ID + "-peaks.tsv")
        or os.path.isfile(experiment_ID + "-narrowpeak.bed")
        or os.path.isfile(experiment_ID + "-summits.bed")
        or os.path.isfile(experiment_ID + "-log.txt")
        or os.path.isfile(experiment_ID + "-model.pdf")
        or os.path.isfile(experiment_ID + "_peaks.narrowpeak")
        or os.path.isfile(experiment_ID + "_peaks.xls")
        or os.path.isfile(experiment_ID + "_summits.bed")
        or os.path.isfile(experiment_ID + "_model.r")
    )
    if exp_bool:
        print(experiment_ID, "exists")
        sys.stdout.flush()
        return True
    return False

def write_bad_exp(exp_accession):

    """
    Logs bad experiments with some odd properties to be revisited later
    """

    f = open("bad_experiments.txt", "a")
    f.write(str(exp_accession))
    f.close()


# Main Code: Requires a files.txt file which can be downloaded from the Encode website after selecting by field in the data Matrix


# DOWNLOAD METADATA AND READ IN AS A DATAFRAME **********************************************************************

files_txt_list = []
with open("files.txt") as f:                                # Open files.txt and pull out link to metadata.
    for line in f:
        files_txt_list.append(line.strip())
metadata = files_txt_list[0]

if os.path.isfile('metadata.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
    download_link(metadata)
    meta_table = pd.read_table('metadata.tsv')
else:
    meta_table = pd.read_table('metadata.tsv')               # If metadata.tsv exists, read as dataframe
print(meta_table)

# ******************************************************************************************************************
# DOWNLOAD HG19 REFERENCE GENOME ***********************************************************************************

global reference_genome
reference_genome = "hg19ind"
ref_download_filename = "chromFa.tar.gz"
if os.path.isfile(ref_download_filename) == False:
    ref_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/" + ref_download_filename
    download_link(ref_url)
    subprocess.run(["gunzip", ref_download_filename])
    subprocess.run(["tar", "xkf", ref_download_filename[0:-3]])
    subprocess.run(["cat", "*.fa", ">>", "reference.fa"])
    subprocess.run(["bowtie2-build", "reference.fa", reference_genome])         #Builds reference genome

# *****************************************************************************************************************
# GENERATE EXPERIMENT TO TRIAL MATCHING **************************************************************************

value_counts = meta_table['Experiment accession'].value_counts() # Finds all unique experiments and number of associated files for each
#print(type(valueCounts), len(valueCounts), "Experiments Conducted")

# Generate mapping of all accession IDs per experiment to each transcription factor
global accession_length
accession_length = 11            # Length of file accession IDs (can change over time?)

experiments_dict = {}             # For each key(Experiment Accession ID), there is a corresponding Experiment Object

for exp, val in value_counts.iteritems():

    experiments_dict[exp] = Experiment(exp, None, [], [])                                  # Add experiment object

    for key, replicate in meta_table.loc[meta_table['Experiment accession'] == exp, 'File accession'].iteritems():
        experiments_dict[exp].noncontrols.append(Replicate(replicate, is_control = False))

    for key, control in meta_table.loc[meta_table['Experiment accession'] == exp, 'Controlled by'].iteritems():

        if type(control) is not str:
            continue

        if len(control) > 2*accession_length:                               # Some controls are strangely organized
            all_controls = control.split(',')
            for indiv_control in all_controls:
                accession_only = indiv_control[-(accession_length+1):-1]
                experiments_dict[exp].controls.append(Replicate(accession_only, is_control = True))
        else:
            accession_only = control[-(accession_length+1):-1]
            control_download = accesion_to_url(accession_only)
            experiments_dict[exp].controls.append(Replicate(accession_only, is_control = True))

print(experiments_dict)
print("number of studies", len(experiments_dict))
print(meta_table["Experiment target"].value_counts())
# *****************************************************************************************************************
# CREATE MACS PEAK FOR EACH EXPERIMENT ****************************************************************************

bad_experiments = []
for exp, experiment in experiments_dict.items():
    if experiment_exists(exp, experiment):
        sys.stdout.write(exp)
        continue
    experiment.controls = list(set(experiment.controls))
    if len(experiment.noncontrols) == 0 or len(experiment.controls) == 0:                                     # Check if replicates or controls are empty
        write_bad_exp(exp)
        continue

    experiment.process_experiment()

sys.stdout.write("FINISHED")


# Local testing example experiment

#expToAccession = {}
#expToAccession['ENCSR000DKX'] = [['ENCFF000RPR', 'ENCFF000RPS'], []]
#expToAccession['ENCSR000DNV'] = [['ENCFF000XBD', 'ENCFF000XBE'], ['ENCFF000XGP', 'ENCFF000XGP']]

#


