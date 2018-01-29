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
# TODO Use Hg38 reference genome

distance_threshold_upstream = -1000                         #Permissable peak distances from transcription start site for inclusion
distance_threshold_downstream = 100
class Experiment:

    def __init__(self, accession_id, transcription_factor, noncontrols, controls):

        self.accession_id = accession_id
        self.transcription_factor = transcription_factor
        self.noncontrols = noncontrols
        self.controls = controls
        self.file_state = None
        self.processed_file = None

    def assess_experiment(self):
        #Checks if .bam/peaks are present and stores this as a property
        print("assess_experiment")
        exp_json = get_request(self.accession_id)                           # Gets json for given experiment and converts to dict
        files = exp_json['files']                                           # Pulls out file list from expreiment json
        files_acc = [(file_string[7:])[:-1] for file_string in files]
        file_props = [find_file_props(acc) for acc in files_acc]
        files = self.determine_file_state(file_props)


    def determine_file_state(self, file_props):
        """
        Determines the current level of processing for the experiment(fastq,bam,peak)
        Also returns a list of accessions for all files of the most-processed file type
        """
        print("determine_file_state")
        def is_fully_processed(self, file_props):
            #first determine biological replicates
            largest_group = []

            for ind_file in file_props:
                if ind_file[3] == None:
                    continue
                if len(ind_file[3]) > len(largest_group):
                    largest_group = ind_file[3]

            target_properties = ["bed narrowPeak", "peaks", largest_group]

            for ind_file in file_props:
                print(ind_file, largest_group)
                if ind_file[1:] == target_properties:
                    print("fully processed")
                    self.processed_file = ind_file[0]
                    return True


        fastq_files = []
        bam_files = []
        peak_files = []

        if is_fully_processed(self, file_props):        #Handles the fully processed case
            print("fully processed")
            self.file_state = "bed"


        # [accession, file_type, output_type, br] <- file preps is a list of these lists

        for single_file in file_props:
            file_type = single_file[1]
            if single_file[1] == "fastq" and len(single_file[3]) == 1:
                fastq_files.append(single_file[0])
            elif single_file[1] == "bam" and single_file[2] == "alignments":
                replicate = Replicate(single_file[0], False)
                bam_files.append(replicate)

        print("BAM FILES:", bam_files)
        print("FASTQ FILES", fastq_files)

        if len(bam_files) > 0:
            print("PROCESSING BAMS")
            self.file_state = "bams"
            self.bam_noncontrol_files = bam_files
            self.bam_control_files = self.find_bam_controls()
            return bam_files

        elif len(fastq_files) > 0:
            self.file_state = "fastq"
            return fastq_files

    def find_bam_controls(self):

        print("find_bam_controls")

        # Finds all control .bam files for the given experiment

        response = get_request(self.accession_id)
        controls = response["possible_controls"]
        controls_acc = [(control_acc[13:])[:-1] for control_acc in controls]
        all_control_bams = []
        all_control_replicates = []

        for control_exp in controls_acc:
            control_bams = self.find_bams_in_controlexps(control_exp)
            all_control_bams = all_control_bams + control_bams

        for single_rep in all_control_bams:
            all_control_replicates.append(Replicate(single_rep, True))

        return all_control_replicates

    def find_bams_in_controlexps(self, control_acc):

        # Given a control experiment accession -> outputs all bam control files
        print("find_bams_in_controlexps")
        exp_json = get_request(self.accession_id)                           # Gets json for given experiment and converts to dict
        files = exp_json['files']                                           # Pulls out file list from expreiment json
        files_acc = [(file_string[7:])[:-1] for file_string in files]
        file_props = [find_file_props(acc) for acc in files_acc]

        target_properties = ["bam", "alignments"]

        control_bams = []

        for ind_file in file_props:
            if ind_file[1:3] == target_properties:
                control_bams.append(ind_file[0])

        return control_bams


    def process_experiment(self):

        # Handles if file has already been processed
        if self.processed_file != None:
            download_link(accession_to_url(self.processed_file,".bed.gz"))
            os.rename(self.processed_file + ".bed.gz", self.accession_id + ".bed.gz")
            subprocess.run(["gunzip", self.accession_id + ".bed.gz"])
            print("downloading completed peaks")
            #sys.stdout.flush()

            return

        if self.file_state == "bams":

            try:

                merged_noncontrols = self.create_merged_bam_no_alignment(False)
                merged_controls = self.create_merged_bam_no_alignment(True)
                self.create_macs_peak(merged_noncontrols, merged_controls)

                print(self.accession_id, merged_noncontrols, merged_controls, "THESE ARE BAMS")
                #sys.stdout.flush()

                return

            except:

                merged_noncontrols = self.create_merged_bam(False)
                merged_controls = self.create_merged_bam(True)
                self.create_macs_peak(merged_noncontrols, merged_controls)

        merged_noncontrols = self.create_merged_bam(False)
        merged_controls = self.create_merged_bam(True)
        self.create_macs_peak(merged_noncontrols, merged_controls)

    def create_merged_bam_no_alignment(self, is_control):

    # Given if we are working with reps or controls, takes the list of files and processes from fastq to merged bam
    # depending on the current level of processing

        if is_control:
            current_reps = self.bam_control_files
            merged_filename = "_merged_controls.bam"
        else:
            current_reps = self.bam_noncontrol_files
            merged_filename = "_merged_noncontrols.bam"

        rep_bam_files = []
        for rep in current_reps:
            bam_file = rep.download_only_bam()
            rep_bam_files.append(bam_file)              # List of (accession + .bam) filenames to be merged

        merged_file = self.accession_id + merged_filename

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

    def create_merged_bam(self, is_control):

        # Given if we are working with reps or controls, takes the list of files and processes from fastq to merged bam
        # depending on the current level of processing

        if is_control:
            current_reps = self.controls
            merged_filename = "_merged_controls.bam"
        else:
            current_reps = self.noncontrols
            merged_filename = "_merged_noncontrols.bam"

        rep_bam_files = []
        for rep in current_reps:
            bam_file = rep.download_to_bam()
            rep_bam_files.append(bam_file)              # List of (accession + .bam) filenames to be merged

        merged_file = self.accession_id + merged_filename

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

    def annotate_experiment(self):

        if self.file_state == None:
            experiment_bed_file = self.accession_id + ".bed"
        else:
            experiment_bed_file = self.accession_id + "_summits.bed"

        annotate_out = self.accession_id +  ".txt"
        annotation_file = open(annotate_out, 'w')
        annotate_process = subprocess.Popen(["annotatePeaks.pl", experiment_bed_file, "hg38"], stderr=None, stdout=annotate_out)
        output = annotate_process.communicate()

    def process_annotation(self):

        annotation_df = pandas.read_csv(self.accession_id + ".txt", sep='\t')
        # Insert filtering criteria below
        distance_df = (annotation_df[(annotation_df["Distance to TSS"] > distance_threshold_upstream) \
                            & (annotation_df["Distance to TSS"] < distance_threshold_downstream) \
                            ])

        all_genes = []
        for index, row in distance_df.iterrows():
            all_genes.append(str(row["Gene Name"]))

        write_list = [self.transcription_factor] + all_genes

        out_file = open(self.accession_id + ".tfgenes", 'w')
        out_file.write("\n".join(write_list))


class Replicate:

    def __init__(self, accession_id, is_control):

        self.accession_id = accession_id
        self.is_control = is_control

    def download_to_bam(self):
        """
        Downloads a fastq for a given accession ID, and processes it into a corresponding .bam file
        """

        subprocess.run(["touch", self.accession_id + ".place"])              # Creates a .place file as a placeholder
        download_link(accession_to_url(self.accession_id, ".fastq.gz"))
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

    def download_only_bam(self):

        download_link(accession_to_url(self.accession_id, ".bam"))
        return self.accession_id + ".bam"

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

def accession_to_url(accession, file_tag):
    """
    Converts a file accession to a relevant download URL
    """

    url = 'https://www.encodeproject.org/files/' + accession + '/@@download/' + accession + file_tag
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
        return False


    for noncontrol in experiment.noncontrols:
        if extension_exists(noncontrol.accession_id):
            return True
    for control in experiment.controls:
        if extension_exists(control.accession_id):
            return True

    exp_bool = (os.path.isfile(experiment_ID + "-peaks.tsv")
        or os.path.isfile(experiment_ID + ".place")
        or os.path.isfile(experiment_ID + ".bed.gz")
        or os.path.isfile(experiment_ID + ".bed")
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
        #sys.stdout.flush()
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

if os.path.isfile('metadata.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
    metadata = files_txt_list[0]
    download_link(metadata)
    meta_table = pd.read_table('metadata.tsv')
else:
    meta_table = pd.read_table('metadata.tsv')               # If metadata.tsv exists, read as dataframe
print(meta_table)

# ******************************************************************************************************************
# DOWNLOAD HG19 REFERENCE GENOME ***********************************************************************************

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
            control_download = accession_to_url(accession_only, "fastq.gz")
            experiments_dict[exp].controls.append(Replicate(accession_only, is_control = True))

print(experiments_dict)
print("number of studies", len(experiments_dict))
print(meta_table["Experiment target"].value_counts())
# *****************************************************************************************************************
# CREATE MACS PEAK FOR EACH EXPERIMENT ****************************************************************************

bad_experiments = []
for exp, experiment in experiments_dict.items():
    print(exp)
    if experiment_exists(exp, experiment):
        sys.stdout.write(exp)
        continue
    experiment.controls = list(set(experiment.controls))
    if len(experiment.noncontrols) == 0 or len(experiment.controls) == 0:                                     # Check if replicates or controls are empty
        write_bad_exp(exp)
        continue
    subprocess.run(["touch", experiment.accession_id + ".place"])
    experiment.assess_experiment()
    experiment.process_experiment()
    experiment.annotate_experiment()

sys.stdout.write("FINISHED")


# Local testing example experiment

#expToAccession = {}
#expToAccession['ENCSR000DKX'] = [['ENCFF000RPR', 'ENCFF000RPS'], []]
#expToAccession['ENCSR000DNV'] = [['ENCFF000XBD', 'ENCFF000XBE'], ['ENCFF000XGP', 'ENCFF000XGP']]

#



