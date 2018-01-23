import wget
import os
import pandas as pd
import numpy as np
import subprocess
import urllib
import sys

# Dependencies:
# python libraries above
# wget
# bowtie2
# samtools
# macs2

# Start with files.txt file

# TO-DO
# DONE: Use samtools quickcheck option to check the file format and headers are correct first
# TODO User subprocess to get the exit status of the commands being run to see if they finished without errors.
# TODO Pipe bowtie alignment output directly to samtools view to avoid .sam file intermediary
# TODO Refactor code!

class Experiment:

    def __init__(self, accession_id, transcription_factor, replicate_IDs, control_IDs):

        self.accession_id = accession_id
        self.transcription_factor = transcription_factor
        self.replicate_IDs = replicate_IDs
        self.control_IDs = control_IDs

    def process_experiment(self):


class Replicate:

    def __init__(self, accession_id, transcription_factor, replicate_IDs, control_IDs):

        self.accession_id = accession_id

    def 




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

def DownloadLink(url):
    """
    Download from link url. If there is an error along the way, tries again with continue invoked
    """

    try:
        print("downloading")
        subprocess.run(["wget", "-c", url])
    except:
        DownloadLink(url)

def DeleteFile(filename):
    """
    # Delete a file if it exists. Otherwise do nothing
    """

    if os.path.isfile(filename) == True:
        os.remove(filename)
        sys.stdout.write(filename, "deleted")
    else:
        return

def AccesionToURL(accession):
    """
    Converts a file accession to a relevant download URL
    """

    URL = 'https://www.encodeproject.org/files/' + accession + '/@@download/' + accession + '.fastq.gz'
    return URL

def DownloadToBam(accession):
    """
    Downloads a fastq for a given accession ID, and processes it into a corresponding .bam file
    """

    subprocess.run(["touch", accession+".place"])              # Creates a .place file as a placeholder
    DownloadLink(AccesionToURL(accession))
    downloadFile = accession + ".fastq.gz"
    print("downloaded file")
    subprocess.run(["gunzip", downloadFile])                       # Unzip the .fastq.gz to a fastq
    print("unzipped")
    # Align .fastq and output .sam; .sam is piped into samtools view to output .bam

    samOut = accession + ".sam"
    bamOut = accession + ".bam"
    #subprocess.Popen("bowtie2 -xpv " + referenceGenome + " -US " + downloadFile[0:-3] + " | samtools view -bS - > " + bamOut, shell = True)          
    bamFile = open(bamOut, 'w')
    bowT = subprocess.Popen(["bowtie2", "-x ", referenceGenome, "-U", downloadFile[0:-3]], bufsize = 500*10**6, stderr=None, stdout=subprocess.PIPE)
    #samT = subprocess.Popen(["cat"], stderr=None, stdin=bowT.stdout, stdout=bamFile)
    samT = subprocess.Popen(["samtools", "view", "-bSu", "-"], stderr=None, stdin=bowT.stdout, stdout=bamFile)
    bowT.stdout.close()
    output = samT.communicate()
    bowT.wait()
    '''
    bamFile = open(bamOut, 'w')
    bowT = subproces.Popen(["bowtie2", "-x", "hg19ind", "-U", "ENCFF796VDT.fastq"], stderr=None, stdout=subprocess.PIPE)
    samT = subprocess.Popen(["cat"], stderr=None, stdin=bowT.stdout, stdout=bamFile)
    bowT.stdout.close()
    output = samT.communicate()
    bowT.wait()
    '''
    print("aligned and output as bam")

    DeleteFile(downloadFile)
    #subprocess.run(["samtools", "view" , "-S", "-b", samOut, "-o", bamOut])                          # Conver .sam to .bam
    
    QuickCheck(bamOut)

    DeleteFile(samOut)

    return bamOut

def CheckIfExperimentExists(experimentID, trials):
    """
    Checks if an experiment has already been processed or is being processed
    """
    for trialType in trials:
        for trial in trialType:
            if os.path.isfile(trial+".fastq.gz"): 
                return True
            if os.path.isfile(trial+".fastq"):
                return True
            if os.path.isfile(trial+".sam"):
                return True
            if os.path.isfile(trial+".bam"):
                return True
            if os.path.isfile(trial+".place"):
                return True

    if os.path.isfile(experimentID + "-peaks.tsv"):
        return True
    if os.path.isfile(experimentID + "-narrowpeak.bed"):
        return True
    if os.path.isfile(experimentID + "-summits.bed"):
        return True
    if os.path.isfile(experimentID + "-log.txt"):
        return True
    if os.path.isfile(experimentID + "-model.pdf"):
        return True
    return False
    

# Main Code: Requires a files.txt file which can be downloaded from the Encode website after selecting by field in the data Matrix


# DOWNLOAD METADATA AND READ IN AS A DATAFRAME **********************************************************************

filesTxt = []
with open("files.txt") as f:                                # Open files.txt and pull out link to metadata.
    for line in f:
        filesTxt.append(line.strip())
metadata = filesTxt[0]

if os.path.isfile('metadata.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
    DownloadLink(metadata)
    metaTable = pd.read_table('metadata.tsv')
else:
    metaTable = pd.read_table('metadata.tsv')               # If metadata.tsv exists, read as dataframe
print(metaTable)

# ******************************************************************************************************************
# DOWNLOAD HG19 REFERENCE GENOME ***********************************************************************************

global referenceGenome
referenceGenome = "hg19ind"
refDownloadFileName = "chromFa.tar.gz"
if os.path.isfile(refDownloadFileName) == False:
    refURL = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/" + refDownloadFileName
    DownloadLink(refURL)
    subprocess.run(["gunzip", refDownloadFileName])
    subprocess.run(["tar", "xkf", refDownloadFileName[0:-3]])
    subprocess.run(["cat", "*.fa", ">>", "reference.fa"])
    subprocess.run(["bowtie2-build", "reference.fa", referenceGenome])

# *****************************************************************************************************************
# GENERATE EXPERIMENT TO TRIAL MATCHING **************************************************************************

valueCounts = metaTable['Experiment accession'].value_counts() # Finds all unique experiments and number of associated files for each
#print(type(valueCounts), len(valueCounts), "Experiments Conducted")

# Generate mapping of all accession IDs per experiment to each transcription factor
global accessionLength 
accessionLength = 11            # Length of file accession IDs (can change over time?)

expToAccession = {}             # For each key(Experiment Accession), there is a list with a sublist(replicates) and another sublist(controls)

for exp, val in valueCounts.iteritems():

    expToAccession[exp] = [[],[]]                                   # Each experiment has replicates and controls
    
    for key, val in metaTable.loc[metaTable['Experiment accession']== exp, 'File accession'].iteritems():
        expToAccession[exp][0].append(val)

    for key, val in metaTable.loc[metaTable['Experiment accession']== exp, 'Controlled by'].iteritems():

        if type(val) is not str:
            continue

        if len(val) > 2*accessionLength:
            allControls = val.split(',')
            for control in allControls:
                accessionOnly = control[-(accessionLength+1):-1]
                expToAccession[exp][1].append(accessionOnly)
        else:
            accessionOnly = val[-(accessionLength+1):-1]
            controlDownload = AccesionToURL(accessionOnly)
            expToAccession[exp][1].append(accessionOnly)

print(expToAccession)                                       
print("number of studies", len(expToAccession))
print( metaTable["Experiment target"].value_counts())
# *****************************************************************************************************************
# CREATE MACS PEAK FOR EACH EXPERIMENT ****************************************************************************

expToAccession = {}
expToAccession['ENCSR000DKX'] = [['ENCFF000RPR', 'ENCFF000RPS'], []]
expToAccession['ENCSR000DNV'] = [['ENCFF000XBD', 'ENCFF000XBE'], ['ENCFF000XGP', 'ENCFF000XGP']]


badExperiments = []
for exp, trials in expToAccession.items():
    if CheckIfExperimentExists(exp, trials) == True:
        sys.stdout.write(exp)
        continue
    trials[1] = list(set(trials[1]))
    if len(trials[0]) == 0 or len(trials[1]) == 0:                                     # Check if replicates or controls are empty
        badExperiments.append(exp)
        continue
    sys.stdout.write(str(trials))
    replicateBamFiles = []                                                             # Downloads all replicates for an experiment,
    for replicateAcc in trials[0]:                                                     # converts them to BAM fil es, and merges them into
        bamFile = DownloadToBam(replicateAcc)                                          # a single mergedRepFile
        replicateBamFiles.append(bamFile)

    mergedRepFile = exp + "mergedReplicates.bam" 

    if len(replicateBamFiles) == 1:
        mergedRepFile = replicateBamFiles[0]
    elif len(replicateBamFiles) > 1: 

        subprocess.run(["samtools", "merge", mergedRepFile]+replicateBamFiles)

    for rBF in replicateBamFiles:                           # Delete unneccessary files(replicates)
        DeleteFile(rBF)                                     # ReplicateAccession.bam
        DeleteFile(rBF[:-4]+".sam")                         # ReplicateAccession.sam
        DeleteFile(rBF[:-4]+".fastq")                       # ReplicateAccession.fastq
        DeleteFile(rBF[:-4]+".place")                       # ReplicateAccession.place  

    controlBamFiles = []                                                            # Does the same as replicates but for controls
    for controlAcc in trials[1]:
        bamFile = DownloadToBam(controlAcc)
        controlBamFiles.append(bamFile)
    mergedControlFile = exp + "mergedControls.bam"
    if len(controlBamFiles) == 1:
        mergedControlFile = controlBamFiles[0]
    elif len(controlBamFiles) > 1:
        try:
            subprocess.run(["samtools", "merge", mergedControlFile]+controlBamFiles)
        except:
            sys.stdout.write("Problem merging control bam files")
            sys.exit()
    for cBF in controlBamFiles:                             # Delete unneccessary files(controls)
        DeleteFile(cBF)                                     # ControlAccession.bam
        DeleteFile(cBF[:-4]+".sam")                         # ControlAccession.sam
        DeleteFile(cBF[:-4]+".fastq")                       # ControlAccession.fastq
        DeleteFile(cBF[:-4]+".place")
    try:
        subprocess.run(["macs2", "callpeak", "-t", mergedRepFile, "-c", mergedControlFile, "-n", exp]) # Call Macs Peak of the merged files
    except:
        sys.stdout.write("Problem running MACS2")
        sys.exit()

    for rBF in replicateBamFiles:                           # Delete unneccessary files(replicates)
        DeleteFile(rBF)                                     # ReplicateAccession.bam
        DeleteFile(rBF[:-4]+".sam")                         # ReplicateAccession.sam
        DeleteFile(rBF[:-4]+".fastq")                       # ReplicateAccession.fastq
        DeleteFile(rBF[:-4]+".place")                       # ReplicateAccession.place
    for cBF in controlBamFiles:                             # Delete unneccessary files(controls)
        DeleteFile(cBF)                                     # ControlAccession.bam
        DeleteFile(cBF[:-4]+".sam")                         # ControlAccession.sam
        DeleteFile(cBF[:-4]+".fastq")                       # ControlAccession.fastq
        DeleteFile(cBF[:-4]+".place")                       # ControlAccession.place
    DeleteFile(mergedRepFile)                               # ExperimentAccessionMergedReplicates.bam
    DeleteFile(mergedControlFile)                           # ExperimentAccessionMergedControls.bam

sys.stdout.write("FINISHED")


# Local testing example experiment

#expToAccession = {}
#expToAccession['ENCSR000DKX'] = [['ENCFF000RPR', 'ENCFF000RPS'], []]
#expToAccession['ENCSR000DNV'] = [['ENCFF000XBD', 'ENCFF000XBE'], ['ENCFF000XGP', 'ENCFF000XGP']]

# 



