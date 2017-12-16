# -*- coding: utf-8 -*-


import sys
import os


from os import listdir
from os.path import isfile, join

import json

import statistics

from pymongo import MongoClient

def filterMinus999(tab):
  newTab = []
  for x in tab:
    if(x != -999):
      newTab.append(x)
  return newTab

def filterDash(tab):
  tabtemp = []
  for x in tab:
    if(x[0] != "-"):
      tabtemp.append(x)
  return tabtemp


print("startLoaderCompile")
client = MongoClient()
client = MongoClient('localhost', 27027)
db = client.rdv
print("connected")

publicFolder = "/u/malricp/rdv/public/"

pathTab = []
pathTab.append({"name":"so_detail","path":["localNcmD_detail","so"],"soft":"so"})
pathTab.append({"name":"mcff_detail","path":["localNcmD_detail","mcff"],"soft":"mcff"})
pathTab.append({"name":"so_ncm","path":["localNcmD","so"],"soft":"so"})
pathTab.append({"name":"mcff_ncm","path":["localNcmD","mcff"],"soft":"mcff"})


di = '/u/malricp/rdv/public/JSON_FOLDER_test/'
dirs = [o for o in os.listdir(di) 
                    if o.startswith("ETERNA")]

for folderName in dirs:
  onlyfiles = [f for f in listdir(di+folderName) if isfile(join(di+folderName, f))]

  #onlyfiles = [file for file in listdir(mypath) if isfile(join(mypath, file))]

  rnaTab = []
  counter = 0

  for file in onlyfiles:
      if(file.endswith(".json") and not file.endswith("file_id_short.json")):
        print("i : "+str(counter))
        counter += 1
        if(counter > 100000):
          break
        rna = json.loads(open(di+folderName+"/"+file,"r").read())
        for dPath in pathTab:
          d = {}
          rnaD = {}
          t = rna["nts"]
          filtered = filterMinus999(rna["scoreTab"])
          nts_score_mean =  statistics.mean(filtered)
          nts_score_sd = statistics.stdev(filtered)
          hi_threshold = nts_score_mean + nts_score_sd
          low_threshold = nts_score_mean 
          for nt in t:
            o = nt
            for p in dPath["path"]:
              #print("keys : "+str(o.keys()))
              o = o[p]
            tab = o.keys()
            for v in tab:
                url = "http://majsrv1.iric.ca:3000/"+folderName+"|"+str(rna["rna_id"])+"|"+str(nt["position"])
                #tbTemp.append(str(v))
                root = v
                freq = str(o[v])
                score = str(nt["score"])
                so_pairee = str(nt["freqPairee_so"])
                mcff_pairee = str(nt["freqPairee_mcff"])
                if(nt["score"] < low_threshold and nt["score"] != -999 and o[v] > 0.1 and o[v] > 0.1):
                  label ="Low"
                  db.ncm_50.insert({"ncm":root,"soft":dPath["soft"],"url":url,"freq":freq,"score":score,"so_pairee":so_pairee,"mcff_pairee":mcff_pairee,"label":label})
                
                if(nt["score"] > hi_threshold):
                  label = "Hi"
                  db.ncm_50.insert({"ncm":root,"soft":dPath["soft"],"url":url,"freq":freq,"score":score,"so_pairee":so_pairee,"mcff_pairee":mcff_pairee,"label":label})
                  
                if((nt["score"] < hi_threshold and nt["score"] != -999) and nt["score"] > low_threshold and o[v] > 0.1 ):
                  label = "Bg"
                  db.ncm_50.insert({"ncm":root,"soft":dPath["soft"],"url":url,"freq":freq,"score":score,"so_pairee":so_pairee,"mcff_pairee":mcff_pairee,"label":label})
            
          
   
  