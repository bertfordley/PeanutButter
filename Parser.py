# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:14:19 2017

@author: rob moseley

Parse through query html files from PeanutButter and create a master table

"""

import os, re
from bs4 import BeautifulSoup
import pandas as pd
from itertools import izip
from joblib import Parallel, delayed
import multiprocessing

def Parser(directory, inputFile, columns):

    cols = columns    
    
    subTable = pd.DataFrame(columns = cols)
    
    # get species name
    species = directory.split("/")[-1]
    
    # create beautifulsoup object of html file
    with open(os.path.join(directory, inputFile), 'r') as htmlFile:
        parsed_html = BeautifulSoup(htmlFile)   
    # get query ID
    queryID = parsed_html.h3.string.split(" ")[3][5:].encode('utf-8')      

    dataMatrix = list()
    
    # get all hits
    ps = parsed_html.find_all('p', {'style' : 'margin-top: 1em; margin-bottom: 0em;'})
    total_ps = len(ps)
    total_hits = parsed_html.find_all('p')[1].text.split(" ")[1]
    if total_hits == "no":
        print("No hits for %s" % queryID)
        total_hits = 0
        total_ps = 0
        percentage = 0
        ID = "NA"
        Hit_Species = "NA"
        Desc = 'NA'
        IDP = 'NA' 
        CP = 'NA'
        Func = "NA"
        Sub = "NA"
        DP = "NA"
        artLink = "NA"
        pubmedID = "NA"
        pubmedLink = "NA"
        More = "NA"
        articleTitle = "NA"
        journal = "NA"
        pubYear = "NA"
        pubmedLink = "NA"
        pubmedID = "NA"
        keyWord = "NA"
        data = [species, queryID.encode('utf-8'), ID.encode('utf-8'), Hit_Species.encode('utf-8'), 
                                           IDP.encode('utf-8'), CP.encode('utf-8'), 
                                           Desc, Func.encode('utf-8'), 
                                           Sub, DP, articleTitle,
                                           journal, pubYear, artLink,
                                           pubmedID.encode('utf-8'), 
                                           pubmedLink.encode('utf-8'), keyWord, More,
                                           total_hits, total_ps, percentage]
        dataMatrix.append(data)
    else:
        percentage = (float(total_ps)/float(total_hits))*100
        percentage = "%.2f" % percentage
        print("Parsing %s: Total Hits = %s; Total p's = %d" % (queryID, total_hits, total_ps))        
    
    
    # keywords = stomatal, stomata, guard cell, guard cells, 
    # keyword regex
    patterns = re.compile(r'\bstomata\b|\bstomatal\b|\bguard cell\b|\bguard cells\b')
    pubmed = re.compile(r'\bPubMed\b')
    # loop through hits and extract info      
    
    for p in ps:       
        
        ID = "NA"
        Hit_Species = "NA"
        Desc = 'NA'
        IDP = 'NA' 
        CP = 'NA'
        Func = "NA"
        Sub = "NA"
        DP = "NA"
        artLink = "NA"
        pubmedID = "NA"
        pubmedLink = "NA"
        More = "NA"
        
        for a in p.find_all('a'):
            if str(a.get('title')) == "None":
                continue
            if a.get('title')[0].isdigit():
                [IDP, CP] = a.string.split(",")
            if a.get('title') in ["SwissProt", "MicrobesOnline", "RefSeq"]:
                ID = a.string.encode('utf-8')
                if a.get('title') == "SwissProt":
                    Desc = p.b.text.encode('utf-8')
                if a.get('title') in ["MicrobesOnline", "RefSeq"]:
                    if ID == "NA":
                        Desc = p.text.partition(" ")[2].encode('utf-8')
                    else:
                        Desc = Desc +  ";" + p.text.partition(" ")[2].encode('utf-8')
                        
        for it in p.find_all('i'):
            Species = it.text.encode('utf-8')
            break
        
        for li in p.next_sibling:
            if not isinstance(li, unicode):
                if li.text == "More":
                    More = li.text

            else:
                articleTitle = ""
                journal = ""
                pubYear = ""
                keyWord = ""
                More = ""
                data = [species, queryID.encode('utf-8'), ID.encode('utf-8'), Hit_Species.encode('utf-8'), 
                                           IDP.encode('utf-8'), CP.encode('utf-8'), 
                                           Desc, Func.encode('utf-8'), 
                                           Sub, DP, articleTitle,
                                           journal, pubYear, artLink,
                                           pubmedID.encode('utf-8'), 
                                           pubmedLink.encode('utf-8'), keyWord, More,
                                           total_hits, total_ps, percentage]
                dataMatrix.append(data)

        for li in p.next_sibling:
            if not isinstance(li, unicode):
                if li.text != "More":  
                    if li.text.startswith("FUNC") or li.text.startswith("SUB") or li.text.startswith("DIS"):                       
                        subDesc = li.text.split("\n")
                        for text in subDesc:
                            if text.startswith("FUNC"):
                                Func = subDesc[0][10:]
                            if text.startswith("SUB"):
                                Sub = text[9:]
                            if text.startswith("DIS"):
                                DP = text[22:]
                    else: 
                        articleTitle = "NA"
                        journal = "NA"
                        pubYear = "NA"
                        pubmedLink = "NA"
                        pubmedID = "NA"
                        keyWord = "NA"
                        for a,i in izip(li.find_all('a'), li.find_all('small')):
                            artLink = a.get("href")
                            if patterns.search(str(li.li)):
                                keyWord = "Stomatal"
                            articleTitle = a.get_text()
                            if i is not None:
                                if pubmed.search(str(i)):
                                    id_link = i.a.next_sibling.next_sibling.get("href")
                                    pubmedLink = id_link
                                    pubmedID = id_link[-8:]
                                    pubmedSplit = i.text.partition(",")[2]
                                    journal = pubmedSplit[1:-14]
                                    pubYear = pubmedSplit[-13:-9]
                                else:
                                    textSplit = i.text.partition(",")[2]
                                    journal = textSplit[1:-6]
                                    pubYear = textSplit[-5:]
                        articleTitle = articleTitle.encode('utf-8')  
                        data = [species, queryID.encode('utf-8'), ID.encode('utf-8'), Hit_Species.encode('utf-8'), 
                                           IDP.encode('utf-8'), CP.encode('utf-8'), 
                                           Desc, Func.encode('utf-8'), 
                                           Sub, DP, articleTitle,
                                           journal, pubYear, artLink,
                                           pubmedID.encode('utf-8'), 
                                           pubmedLink.encode('utf-8'), keyWord, More,
                                           total_hits, total_ps, percentage]
                        dataMatrix.append(data)
            
    subTable = pd.DataFrame(dataMatrix, columns=cols)
    
    outfile = directory + "_output.txt"
    
    if os.stat(outfile).st_size == 0:
        subTable.to_csv(outfile, mode="a", header=True, sep="\t", index=False, encoding="utf-8")  
    else:
        subTable.to_csv(outfile, mode="a", header=False, sep="\t", index=False, encoding="utf-8")

if __name__ == "__main__":
    
    n_jobs = multiprocessing.cpu_count()    
    
    cols = ["Species",  "Query_ID",  "ID", "Hit_Species", "Identity_%",  "Coverage_%", 
    "Description", "Function", "Subunit", "Disruption_Phenotype", 
    "Article_Title", "Journal", "Publication_Year", "Journal_Link", 
    "PubMed_ID", "PubMed_Link", "Keyword_Match", "More_Link", "Total_Hits",
    "Total_Hits_Listed", "Percent_Covered"]
      
    topDir = "/Volumes/RMdata1A/PeanutButter20170711/output_html/media/y9x/XY8T02A/vip/PeanutButter_search/out" 

# Non-parallel
#    for spDir in os.listdir(topDir):
##        if not spDir.startswith('.') and (spDir == "Kafe" or spDir == "Arth"):
#        if not spDir.startswith('.') and (spDir == "test1" or spDir == "test2"): 
#            OUT = open(os.path.join(topDir, spDir + "_output.txt"), "w")
#            OUT.close()
#            head = False
#            directory = os.path.join(topDir, spDir)
#            for htmlfile in os.listdir(directory):
#                if not htmlfile.startswith("."):
#                    print(htmlfile)
#                    species_df = Parser(directory, htmlfile, cols)
#                    if not head:
#                        species_df.to_csv(os.path.join(topDir, spDir + "_output.txt"),
#                                      mode="a", header=True, sep="\t", index=False)
#                        head = True
#                    else:
#                        species_df.to_csv(os.path.join(topDir, spDir + "_output.txt"),
#                                      mode="a", header=False, sep="\t", index=False)

    for spDir in os.listdir(topDir):
#        if not spDir.startswith('.') and (spDir == "Kafe" or spDir == "Arth"):
        if not spDir.startswith('.') and (spDir == "test3"):# or spDir == "test2"): 
            OUT = open(os.path.join(topDir, spDir + "_output.txt"), "w")
            OUT.close()
            head = False
            directory = os.path.join(topDir, spDir)
            Parallel(n_jobs=n_jobs) (delayed(Parser) (directory, htmlfile, cols) for htmlfile in os.listdir(directory))
