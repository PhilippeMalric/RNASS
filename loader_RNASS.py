#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import RNASS
import json
import math
import handler
from os import listdir
from os.path import isfile, join



def scoreArrayAug(scoreA,index,f):
  if(len(scoreA) == 1):
    return [0]
  for i in range(0,int(index)-1):
    scoreA.insert(0,-999)
  while(len(scoreA) < f):
    scoreA.append(-999)  
  if(scoreA[-1] == 0):
    score = 0
    i = -1
    while(score == 0):
      scoreA[i] = -999
      i = i-1
      if(len(scoreA)>1 and i*-1 < len(scoreA)):
        score = scoreA[i]
      else: 
        return scoreA
  return scoreA[0:f]


'''
f = open("/u/malricp/_MCfoldVsRNAsubopt/projet/Mam4/folder2.txt",'r')
r = f.read()
f.close()
folders = r.split("\n")
'''


id_num = 0
publicFolder = "/u/malricp/rdv/public/"

'''
lineId = int(sys.argv[1])
lineId2 = int(sys.argv[2])
file = sys.argv[3]
'''

lineId = 0
lineId2 = 1
file = "/u/malricp/_MCfoldVsRNAsubopt/projet/Mam4/rmdb/ETERNA_R82_0000.rdat"


stnSeuil = 1
stn = 0

rna2DTab = []

folder_in = "rmdb"
linesParsed = []
try:
  root = file.split("/")[-1][:-5]
  rdat = handler.RDATFile()
  rdat.load(open(file))
  rdat.validate()
  global_annotation = rdat._annotation_str(rdat.annotations,"|")
  for (key,value) in rdat.constructs.items() :
    if(hasattr(rdat.constructs[key].data[0],"annotations")):
      print("annotation in data")
      print("annotations : "+str(rdat.constructs[key].data[0].annotations))
      annot_tab_str = [x.annotations for x in rdat.constructs[key].data]
      if("sequence" in rdat.constructs[key].data[0].annotations):
          print ("sequence in annotation")
          rdat.seqTab = [x.annotations["sequence"][0] for x in rdat.constructs[key].data]
          #print("seqTab : "+"\n".join(rdat.seqTab))
          print("seqTab len : "+str(len(rdat.seqTab)))
          print("lineId  :"+str(lineId) + " len(rdat.seqTab) : "+str(len(rdat.seqTab)) + " len(value.data) : "+str(len(value.data)))
      if(lineId < len(rdat.seqTab) and lineId < len(value.data) ):
        if(lineId2 > len(value.data)):
          lineId2 = len(value.data)
        for i in range(lineId,lineId2):
          seq = rdat.seqTab[i]
          v = value.data[i].values
          e = value.data[i].errors
          offS = rdat.constructs[key].offset
          v = scoreArrayAug(v,offS,len(seq))
          e = scoreArrayAug(e,offS,len(seq))
          
          if("signal_to_noise" in rdat.constructs[key].data[0].annotations):
            print ("signal_to_noise in annotation")
            stn = rdat.constructs[key].data[i].annotations["signal_to_noise"][0].split(":")[1]
          else : 
            print ("No signal_to_noise in annotation")
            stn = 0
          '''
          print("i :"+str(i)+" seq : "+seq)
          print("offset : "+str(rdat.constructs[key].offset))
          print("len seq : "+str(len(seq)) +" v : "+str(len(v)))
          '''
          
          oneLineParsed =  {}
          oneLineParsed.clear()
          oneLineParsed["stn"] = stn
          oneLineParsed["seq"] = seq
          oneLineParsed["structure"] = ""
          oneLineParsed["scoreTab"] =  v
          oneLineParsed["errorTab"] = e
          oneLineParsed["name"] = ""
          oneLineParsed["title"] = ""
          oneLineParsed["type_id"] = ""
          oneLineParsed["nid"] = "1"
          oneLineParsed["projetTitle"] = ""
          oneLineParsed["date"] = ""
          oneLineParsed["autor"] = ""
          oneLineParsed["i"] = i
          oneLineParsed["local_annotation"] = rdat.constructs[key].data[i].annotations
          oneLineParsed["global_annotation"] = global_annotation
          print("stn : "+str(stn))
          if (float(stn) > stnSeuil):
             linesParsed.append(oneLineParsed)
except ValueError:
  print ("erreur!")
  print(ValueError)



print("linesParsed : "+str(len(linesParsed)))
          
if(len(linesParsed)>0):
  id_num += 1
  filePath  = "/u/malricp/_MCfoldVsRNAsubopt/projet/Mam4/"
  #RNA_featureChoice.RNA2D(filePath,root,linesParsed,folder,"experience_6_nov_2017",id_num) 
  d = {
      'inputFileLoader': root,
      'filePath': filePath,
      'folder_in': folder_in,
      'publicFolder': publicFolder,
      'currentExp': "experience_15_dec_2017",
      'uniqueId': lineId,
      'soTreshold': 10,
      'mcffTreshold':10,
      'RNAJSON_in':"JSON_FOLDER_test/"+root,
      'RNAJSON':"JSON_FOLDER_test/"+root,
      'CSVFolder':"CSV_FOLDER/"+root,
      'prediction': False,
      'printCsv':False,
      'so_e_value':10,
      'lowScore' : 0.5,
      'hiScore' : 1,
      'predCutOff' : 0.5,
      'stnSeuil' : stnSeuil,
      'nbNt_Ncm' : 10  
      }
  
  
  #RNASS_rdat_qsub.RNASS(d,[linesParsed[-1]])
  RNASS.RNASS(d,linesParsed)

    
    
      #line = s+";"+v+";"+e+";"+name+";"+title+";"+type_id+";"+nid+";"+projetTitle +";"+ date
                            
 
'''
outPut = "["+",".join([x.basicD for x in rna2DTab])+"]"
f = open("/u/malricp/rdv/public/js/allSmall2.json",'w')
f.write(outPut)
f.close()
'''
