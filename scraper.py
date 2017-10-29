import requests
from bs4 import BeautifulSoup
import wget
import os

# Fand all URLs from website

def PullLinks(url):

    response = requests.get(url)
    # parse html
    page = str(BeautifulSoup(response.content))
    linksList = []

    while True:
        url, n = GetURL(page)
        page = page[n:]
        if url:
            linksList.append(url)
        else:
            break

    return linksList

def GetURL(page):
    """

    :param page: html of web page (here: Python home page) 
    :return: urls in that page 
    """
    start_link = page.find("a href")
    if start_link == -1:
        return None, 0
    start_quote = page.find('"', start_link)
    end_quote = page.find('"', start_quote + 1)
    url = page[start_quote + 1: end_quote]
    return url, end_quote

# Filters 

def FileTypeFilter(links, filetype):

    filteredLinks = []

    for labInd in range(len(links)):
        filteredLabLinks = []
        for labLink in links[labInd]:
            if len(labLink) > len(filetype):
                if labLink[-len(filetype):] == filetype:
                    filteredLabLinks.append(labLink)
        filteredLinks.append(filteredLabLinks)
    
    return filteredLinks

def DownloadLink(url):

    currentFile = wget.download(url)
    print(currentFile, "downloaded")
    return currentFile

def DeleteFile(filename):

    os.remove(filename)
    print(filename, "deleted")

def FileTxtToDict(mainURL, metadataList, dictionary, filetype):
    """
    Reads in each files.txt file. Then adds the metadata for each .fastq.gz file to the main storage dictionary
    """
    DownloadLink(mainURL + "files.txt")

    with open("files.txt") as f:

        for line in f:
            
            tabSplit = line.split("\t")
            scolonSplit = tabSplit[1].split("; ")
            cleanMeta = ['']*len(metadataList)

            for element in scolonSplit:
                equalSplit = element.split("=")

                for i in range(len(metadataList)):
                    if equalSplit[0] == metadataList[i]:
                        cleanMeta[i] = equalSplit[1]

            if tabSplit[0][-len(filetype):] == filetype:
                dictionary[tabSplit[0]] = cleanMeta

    DeleteFile("files.txt")

# List of relevant URLs to be scraped

filetype = ".fastq.gz"                                      # Filetype to be scraped

allUrls = ["http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUchicagoTfbs/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromChip/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/"]

linksForUrls = []       # storage of links 

# Initialize dictionary "table" with metadata headers
metadataList = ["lab", "cell", "antibody", "dateSubmitted", "dateUnrestricted", "peaks"]
dataTable = {}
dataTable["Headers"] = metadataList

# Collects all links from each webpage, and then filters them for only the filetypes of interest

for url in allUrls:

    FileTxtToDict(url, metadataList, dataTable, filetype)         # 

    links = PullLinks(url)
    linksForUrls.append(links)

filteredLinks = FileTypeFilter(linksForUrls, filetype)

print(dataTable)

# Loops through each file link and processes them

# for labUrlInd in range(len(filteredLinks)):
#     for fileLink in filteredLinks[labUrlInd]:

#         newRow = ["", "", "", ]

#         fullUrl = allUrls[labUrlInd] + fileLink

#         currentFile = DownloadLink(fullUrl)                 # Downloads a given link

#         DeleteFile(currentFile)                             # Deletes file once it has been processed

