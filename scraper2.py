import wget
import os
import pandas as pd
import numpy as np
import subprocess

# Start with files.txt file

def DownloadLink(url):

    currentFile = wget.download(url)
    print(currentFile, "downloaded")
    return currentFile

def DeleteFile(filename):

    os.remove(filename)
    print(filename, "deleted")

def AccesionToURL(accession):

    URL = 'https://www.encodeproject.org/files/' + accession + '/@@download/' + accession + '.fastq.gz'
    return URL

def SingleExperiment(experiment, experimentFiles):

    replicateFiles = []
    controlFiles = []

def DownloadToBam(accession):

    downloadFile = DownloadLink(AccesionToURL(accession))
    subprocess.call(["gunzip", downloadFile])
    samOut = accession + ".sam"
    subprocess.call(["bowtie2", "-x", referenceGenome, "-U", downloadFile[0:-3], "-S", samOut])
    bamOut = accession + ".bam"
    subprocess.call(["samtools", "view" , "-S", "-b", samOut, ">", bamOut])

    return bamOut


# Main Code: Requires a files.txt file which can be downloaded from the Encode website after selecting by field in the data Matrix


# DOWNLOAD METADATA AND READ IN AS A DATAFRAME **********************************************************************

filesTxt = []
with open("files.txt") as f:                                # Open files.txt and pull out link to metadata.
    for line in f:
        filesTxt.append(line.strip())
metadata = filesTxt[0]

if os.path.isfile('metadata.tsv') == False:                 # If metadata.tsv does not exist download it and read it as a dataframe
    metaFile = DownloadLink(metadata)
    metaTable = pd.read_table(metaFile)
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
    refDown = wget(refURL)
    subprocess.call(["tar", "xkf"].append(refDown[0:-3]))
    subprocess.call(["cat", "*.fa", ">>", "reference.fa"])
    subprocess.call(["bowtie2-build", "reference.fa", referenceGenome])

# *****************************************************************************************************************
# GENERATE EXPERIMENT TO TRIAL  MATCHING **************************************************************************

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

# *****************************************************************************************************************
# CREATE MACS PEAK FOR EACH EXPERIMENT ****************************************************************************

for exp, trials in expToAccession.items():
    
    replicateBamFiles = []                                                          # Downloads all replicates for an experiment,
    for replicateAcc in trials[0]:                                                     # converts them to BAM files, and merges them into
        bamFile = DownloadToBam(replicateAcc)                                          # a single mergedRepFile
        replicateBamFiles.append(bamFile)
    mergedRepFile = "mergedReplicates.bam"
    subprocess.call(["samtools", "merge", mergedRepFile]+replicateBamFiles)
    for rBF in replicateBamFiles:
        DeleteFile(rBF)

    controlBamFiles = []                                                            # Does the same as replicates but for controls
    for controlAcc in trials[1]:
        bamFile = DownloadToBam(controlAcc)
        controlBamFiles.append(bamFile)
    mergedControlFile = "mergedControls.bam"
    subprocess.call(["samtools", "merge", mergedControlFile]+controlBamFiles)
    for cBF in controlBamFiles:
        DeleteFile(cBF)

    subprocess.call(["macs2", "callpeak", "-t", mergedRepFile, "-c", mergedControlFile, "-n", exp])


# NEED TO CONVERT TO PYTHON 2.7 ********************************************************






