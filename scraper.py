import requests
from bs4 import BeautifulSoup
import wget
import os

# FINDS ALL URLS FOR WEBSITE


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

def FileTypeFilter(links, filetype):

    filteredLinks = []

    for labInd in range(len(links)):
        filteredLabLinks = []
        for labLink in links[labInd]:
            if len(labLink) > len(filetype:
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


# List of relevant URLs to be scraped

allUrls = ["http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUchicagoTfbs/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromChip/",
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwTfbs/"]

linksForUrls = []

# Collects all links from each webpage, and then filters them for only the filetypes of interest

for url in allUrls:
    linksForUrls.append(PullLinks(url))

filteredLinks = FileTypeFilter(linksForUrls, ".fastq.gz")
print(filteredLinks)

for labUrlInd in range(len(filteredLinks)):
    for fileLink in filteredLinks[labUrlInd]:

        fullUrl = allUrls[labUrlInd] + fileLink

        currentFile = DownloadLink(fullUrl)

        DeleteFile(currentFile)