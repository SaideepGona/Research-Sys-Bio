import requests
from bs4 import BeautifulSoup
import wget
import os
import pandas as pd
import numpy as np

# Start with files.txt file

def DownloadLink(url):

    currentFile = wget.download(url)
    print(currentFile, "downloaded")
    return currentFile

def DeleteFile(filename):

    os.remove(filename)
    print(filename, "deleted")


# Main

filesTxt = []
with open("files.txt") as f:
    for line in f:
        filesTxt.append(line.strip())
print(filesTxt)
metadata = filesTxt[0]

metaFile = DownloadLink(metadata)

metaTable = pd.read_table(metaFile)
print(metaTable)

valueCounts = metaTable['Experiment target'].value_counts()
print(type(valueCounts), len(valueCounts), "Unique Transcription Factors Profiled")

for tf, val in valueCounts.iteritems():
    
    for key, val in 


