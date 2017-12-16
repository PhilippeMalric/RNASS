# -*- coding: utf-8 -*-
import statistics
import operator
import os
import subprocess
import itertools
import re
import math
import json
import operator
import urllib
import time
import string
import itertools
import csv
from collections import Counter
from itertools import groupby
from subprocess import Popen,PIPE
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from pymongo import MongoClient

#from scipy import stats

#from Graph import Graph

#from Score import Score


#ecrit des fichier csv avec des choix de caracteristique different
# -2 - position relative
# -1 - sorte seulement
# 0 - subseq only
# 0.1 - pairee ou non
# 1 - ncm only
# 2 - in helice only
# 3 - reduce set of ncmCaracteristic_d
# 4 - All in


class RNASS:
      'generate an html file with d3.js in it to visualize reactivity data on secondary structure'
      count = 0 
      

#filePath,name,seqProfilTab,folder,publicF,expName,id_num
      
      def __init__(self,d,linesParsed):
          
          
          #populate self variable----------------------------------------------
          self.defVar(d,linesParsed)
          
          
          
          self.outputFolder = "js_eternaTsv_new2_number"
          self.outputFolder1 = "js_eternaTsv_new2_number_andPred"
          self.loadCaracteristic()
          
          dataSetTab = []
            
          dataSetTab.append(self.dictForSorte())
          dataSetTab.append(self.dictForBasic())
          dataSetTab.append(self.dictForPairPartnerPathGen())
          dataSetTab.append(self.dictForRegion())
          dataSetTab.append(self.dictForSeqMotif3())
          dataSetTab.append(self.dictForSeqMotif5())
          dataSetTab.append(self.dictForNcm25_mcff())
          dataSetTab.append(self.dictForNcm25_so())
          dataSetTab.append(self.dictForMfeNcm_mcff())
          dataSetTab.append(self.dictForMfeNcm_so())
          dataSetTab.append(self.dictForNcm_mfe_detail_mcff())
          dataSetTab.append(self.dictForNcm_mfe_detail_so())
          dataSetTab.append(self.dictForNcm_detail_mcff())
          dataSetTab.append(self.dictForNcm_detail_so())
          dataSetTab.append(self.dictForNcmVoisin())
          dataSetTab.append(self.dictForHelice())
          
          self.eraseFileAndFillHeader(dataSetTab)
          
          print("start GettingScore")
          self.db = self.openMongoClient()
          print("end GettingScore")
          
          
          if(True or os.path.isfile(self.publicFolder+self.rnaJSONFolder_in+"/"+str(id_num)+".json")):
            
            
            print ("begining of RNA.py -----------------------------")
            RNASS.count += 1
            
            # create a folder where the json files and the csv file will be output
            
            InH_sens = 0.75
            
            
            path = "/u/malricp/rdv/public/"+self.outputFolder+"/"
            if not os.path.exists(path):
                            os.makedirs(path)
            
            #print("self.seqProfilTab : "+str(self.seqProfilTab))
            meanLenSeq = statistics.mean(len(x["seq"]) for x in self.seqProfilTab)
            print("meanLenSeq : "+ str(meanLenSeq))
            
            
            rna_id_num = -1
            nt_id_num = 0
            
            #pour reduire les dimensions
            self.ncmCaracteristicsHi = self.createCaracteristicHi()
            self.ncmCaracteristicsLow = self.createCaracteristicLow()
            
            
            
            
                
            self.cpjp = {}
            
            #self.ncmCaracteristic_d = self.ncmCaracteristicsInit()
            
            #print("self.ncmCaracteristic_d : "+str(len(self.ncmCaracteristic_d.keys())))
            
            fn = "/u/malricp/rdv/public/"+self.outputFolder+"/"+self.name+".csv"
            os.remove(fn) if os.path.exists(fn) else None
            
            '''
            os.remove("./ncm.txt") if os.path.exists("./ncm.txt") else None
            self.fncm = open("./ncm.txt",'w')
            '''
            print("self.seqProfilTab len : " + str(len(self.seqProfilTab)))
            
            for e in self.seqProfilTab:
                    
                    self.seq = e["seq"]
                    self.ms_id =int(self.ms_id)+ 1
                    filtered = self.filterMinus999(e["scoreTab"])
                    if(len(filtered)>1): 
                      self.nts_score_mean =  statistics.mean(filtered)
                      self.nts_score_sd = statistics.stdev(filtered)
                    else:
                      print("scoreTab : "+str(e["scoreTab"]))
                      break
                    if(float(e["stn"]) < self.d['stnSeuil'] or os.path.exists(self.publicFolder+self.rnaJSONFolder+"/"+str(self.ms_id)+".json") or self.nts_score_mean > 1.5 or self.calculFreqR("A",e["seq"]) > 0.5):
                      print("freqA : "+str(self.calculFreqR("A",e["seq"])))
                      print("self.nts_score_mean : "+str(self.nts_score_mean))
                      break
                    
                    lowScore = self.d['lowScore']
                    hiScore =  self.d['hiScore']
                    #print("e['scoreTab'] : "+str(e["scoreTab"]))   
                    self.reactivityVector = []
                    for sc in e["scoreTab"]:
                        #print("sc : "+str(sc))
                        #print("lowScore : "+str(lowScore))
                        if(sc > hiScore):
                            self.reactivityVector.append("Hi")
                        else: 
                            if(sc < lowScore and sc != -999):
                                self.reactivityVector.append("Low")
                            else:
                                if(sc != -999):
                                    self.reactivityVector.append("Bg")
                                else:
                                    self.reactivityVector.append("NA")
                                    
                    #print("self.reactivityVector : "+str(self.reactivityVector))                  
                    print("seqUnderTreatment : "+e["seq"])
                    #print("structure : "+e["structure"])
                    rna_id_num += 1
                    print("rna_id_num : "+str(rna_id_num))
                    
                    #print('e["scoreTab"] = '+str( e["scoreTab"]))#[:40]
                    #print('e["errorTab"] = '+str(e["errorTab"]))#[:40]
                    
                    print('ms_id = '+str(self.ms_id))#[:40]
                    sE_Tab_So = self.rnaSO(e["seq"],True)
                    self.sESTabsubOpt = sE_Tab_So
                    self.tabSubOpt_allRNA_so.append(sE_Tab_So)
                    print("seq : "+e["seq"])
                    sE_Tab_Mcff = self.mcff(e["seq"],True) 
                    print("sE_Tab_Mcff : "+ str(len(sE_Tab_Mcff)))
                    self.sESTabMcff = sE_Tab_Mcff
                    self.tabSubOpt_allRNA_mcff.append(sE_Tab_Mcff)
                    
                    (self.ed,self.freqMfe) = (0,0)
                    
                    #Donne la frequence de pairage des nt
                    frP_mcff = self.sumNtWise([x.struct for x in self.sESTabMcff]) 
                    frP_so = self.sumNtWise([x.struct for x in sE_Tab_So])
                    
                    
                    if(sE_Tab_So and sE_Tab_Mcff):
                        #---------------------------------------------affichage
                        '''
                        print("sE_Tab_SubOpt")
                        print("len sE_Tab_So :" + str(len(sE_Tab_So)))
                        #for s in self.sE_Tab_SubOpt:
                                #print(s.toString())
                        print("sE_Tab_Mcff")
                        print("len sE_Tab_Mcff :" + str(len(sE_Tab_Mcff)))
                        #for s in self.sE_Tab_Mcff:
                                #print(s.toString())
                        '''
                        #-----------------------------------------------------fin affichage
                        
                        
                        #self.nhn_mcff = moyenne(nombre d'helice / longeur)
                        #self.nPbN_mcff = moyenne(nombre de paire de base / longeur)
                        #self.nPbSn_mcff = moyenne(nombre de paire de base stable (ayant toujours le meme partenaire) / longeur)
                        
                        
                        #self.nhn_so = moyenne(nombre d'helice / longeur)
                        #self.nPbN_so = moyenne(nombre de paire de base / longeur)
                        #self.nPbSn_so = moyenne(nombre de paire de base stable (ayant toujours le meme partenaire) / longeur)
                        
                        
                        
                        self.pairTabAllSubOpt_so = [self.findPairGen(x.struct) for x in sE_Tab_So]
                        jsNode_so = [self.jsNodeAll(pairTab,len(e["scoreTab"])) for pairTab in self.pairTabAllSubOpt_so]
                        #print("jsNodePrint_so: "+jsNodePrint_so)
                        
                        self.pairTabAllSubOpt_mcff = [self.findPairGen(x.struct) for x in sE_Tab_Mcff]
                        jsNode_mcff = [self.jsNodeAll(pairTab,len(e["scoreTab"])) for pairTab in self.pairTabAllSubOpt_mcff]
                        #print("jsNodePrint_mcff: "+jsNodePrint_mcff)
                        
                        
                        (t1_mcff,mfp_mcff) = self.findMostfreqPartner(self.pairTabAllSubOpt_mcff )
                        (t1_so,mfp_so) = self.findMostfreqPartner(self.pairTabAllSubOpt_so)
                        
                        #--------------------------------------------------------------
                        (self.ncmD,self.mcffNts,self.soNts,self.colorPredTab_so,self.colorPredTab_mcff,ncmTabDG_so,ncmTabDG_mcff,vide,vide) = self.pbToNcmAll(len(e["seq"]),False)
                        (self.ncmD_mfe,self.mcffNts_mfe,self.soNts_mfe) = self.pbToNcmAll_mfe(len(e["seq"]),False)
                        
                        (self.ncmD_detail,self.mcffNts_detail,self.soNts_detail,self.colorPredTab_so,self.colorPredTab_mcff,ncmTabDG_so,ncmTabDG_mcff,sc_so,sc_mcff) = self.pbToNcmAll(len(e["seq"]),True)
                        (self.ncmD_mfe_detail,self.mcffNts_mfe_detail,self.soNts_mfe_detail ) = self.pbToNcmAll_mfe(len(e["seq"]),True)
                        #-----------------------------------------------------------
                        
                        sc_mcff_t = self.transposeSc(sc_mcff)
                        sc_so_t = self.transposeSc(sc_so)
                        
                        #print("tabVoisinAllsub_so[10] : "+str(self.tabVoisinAllsub_so[10]))
                        #print("tabVoisinAllsub_mcff[10] : "+str(self.tabVoisinAllsub_mcff[10]))
                        
                        dotBraquetMcffTab = [x.struct for x in sE_Tab_Mcff]
                        dotBraquetSoTab = [x.struct for x in sE_Tab_So ]
                        
                      
                        rna_id_num_str = str(rna_id_num)
                        seq = e["seq"]
                        
                        
                        dotBr2p = self.jsDictGeneratorByPredictor(dotBraquetMcffTab,dotBraquetSoTab)
                        d3ForceLayout2p = self.jsDictGeneratorByPredictor(jsNode_mcff,jsNode_so)
                        
                        
                        _2NtMotif = self.extractMotif(2,seq)
                        #_2NtMoti#file = "Mapseek1.txt"f_str = json.dumps(_2NtMotif)
                        scoreTabStr = "["+",".join(['%.5f' % x for x in e["scoreTab"]])+"]"
                        errorTabStr = "["+",".join(['%.5f' % x  for x in e["errorTab"]])+"]"

                        nTs =  []
                        # pour tous les nucleotide -------------------------------
                        for i in range(0,len(seq)):
                            voisinPairAllsub_so = [x[i] for x in self.tabVoisinAllsub_so]
                            voisinPairAllsub_mcff = [x[i] for x in self.tabVoisinAllsub_mcff]
                            #print("voisinPairAllsub_so[0]  : "+str(voisinPairAllsub_so[0]))
                            #print("voisinPairAllsub_mcff[0] : "+str(voisinPairAllsub_mcff[0]))
                            nt_id_num = nt_id_num + 1
                            score = e["scoreTab"][i]
                            erreur = e["errorTab"][i]
                            sorte = e["seq"][i]
                            
                            if(not(i< 1 or i > len(seq)-2)):
                              seqMotif3 = e["seq"][i-1:i+2]
                            else:
                              seqMotif3 = "-"
                            if(not(i< 2 or i > len(seq)-3)):
                              seqMotif5 = e["seq"][i-2:i+3]
                            else:
                              seqMotif5 =  "-"
                            
                            freqPairee_mcff = frP_mcff[i]
                            freqPairee_so = frP_so[i]
                            position = i
                            id_str = str(nt_id_num)
                            etatMcff = [x.struct[i] for x in sE_Tab_Mcff]
                            etatSo = [x.struct[i] for x in sE_Tab_So]
                            etat = self.jsDictGeneratorByPredictor(sum([self.convertDBTo01(x) for x in etatMcff])/float(len(etatMcff)),sum([self.convertDBTo01(x) for x in etatSo])/float(len(etatSo)))
                            self.localNcmD = self.localNcmD_gen(self.mcffNts,self.soNts,i)
                            self.localNcmD_mfe = self.localNcmD_gen(self.mcffNts_mfe,self.soNts_mfe,i)
                            self.localNcmD_detail = self.localNcmD_gen(self.mcffNts_detail,self.soNts_detail,i)
                            self.localNcmD_mfe_detail = self.localNcmD_gen(self.mcffNts_mfe_detail,self.soNts_mfe_detail,i)
                            reactivity_pred = self.reactivityVector[i]
                            #
                            ncmTabDG_so_i = ncmTabDG_so[i]
                            ncmTabDG_mcff_i = ncmTabDG_mcff[i]
                            inHeliceTemp = 0
                            if(i>0):
                              
                                #savoir si le nucleotide est au milieu d'une suite de 3 nucleotide pairee
                                if("2_2" in nTs[i-2]["localNcmD"]["so"] and "2_2" in self.localNcmD["so"] ):
                                  #print("2_2_so both side")
                                  if(float(nTs[i-2]["localNcmD"]["so"]["2_2"]) > InH_sens and  float(self.localNcmD["so"]["2_2"])> InH_sens ):
                                    #print("Hier than :"+ str(InH_sens))
                                    inHeliceFinal = 1
                                  else:
                                    inHeliceFinal = 0
                                else:
                                  #print("keys : "+str(nTs[i-2]["localNcmD"]["so"].keys()))
                                  inHeliceFinal = 0
                              
                              
                                voisinG = self.createVoisin(nTs[i-1])
                                nt_i = i
                                l = len(seq)
                                region = self.findRegion(nt_i,l)
                                nt = self.createNt(id_str,voisinG,erreur,score,sorte,position,etat,self.localNcmD,self.localNcmD_mfe,self.localNcmD_detail,self.localNcmD_mfe_detail,self.ms_id,inHeliceTemp,freqPairee_mcff,freqPairee_so,seqMotif3,seqMotif5,mfp_mcff[i],mfp_so[i],region,voisinPairAllsub_so,voisinPairAllsub_mcff,self.colorPredTab_so[i],self.colorPredTab_mcff[i],reactivity_pred,ncmTabDG_so_i,ncmTabDG_mcff_i,sc_mcff_t[i],sc_so_t[i],self.d["inputFileLoader"])
                                voisinD = self.createVoisin(nt)
                                nTs[i-1]["voisinDroit"] = voisinD
                                
                                
                                  
                                nTs[i-1]["inHelice"] = inHeliceFinal
                            else:
                                voisinG = None
                                inHeliceTemp = 0
                            
                            l = len(seq)
                            region = self.findRegion(i,l)    
                            nt = self.createNt(id_str,voisinG,erreur,score,sorte,position,etat,self.localNcmD,self.localNcmD_mfe,self.localNcmD_detail,self.localNcmD_mfe_detail,self.ms_id,inHeliceTemp,freqPairee_mcff,freqPairee_so,seqMotif3,seqMotif5,mfp_mcff[i],mfp_so[i],region,voisinPairAllsub_so,voisinPairAllsub_mcff,self.colorPredTab_so[i],self.colorPredTab_mcff[i],reactivity_pred,ncmTabDG_so_i,ncmTabDG_mcff_i,sc_mcff_t[i],sc_so_t[i],self.d["inputFileLoader"])
                            nTs.append(nt)
                            if (i == len(seq)-1):
                                nTs[i]["voisinDroit"] = None
                            #if(i == self.testNCM):
                              #print("localNcmD : "+str(localNcmD))
                            #if(nt["score"] != -999):
                              #self.db.ntP_10_Poids.insert({"id":str(self.d["inputFileLoader"])+"|"+str(self.ms_id)+"|"+str(i),"freqPairee_mcff":freqPairee_mcff,"freqPairee_so":freqPairee_so,"label":reactivity_pred})
                              #self.db.nt_10_fin2.insert(nt.copy())
                              
                              
                        #fin pour tous les nucleotide ----------------------------
                        scoreT = [x["score"] for x in nTs]
                        #try:
                        filtered = self.filterMinus999(scoreT)
                        if(len(filtered)>1): 
                          self.nts_score_mean =  statistics.mean(filtered)
                          self.nts_score_sd = statistics.stdev(filtered)
                        else:
                          print("score T  : "+str(scoreT))
                          break
                        lowScore = self.d["lowScore"]
                        hiScore = self.d["hiScore"]
                            
                        nTs_Hi = []
                        nTs_Low = []
                        nTs_Bg = []
                        for nt in nTs:
                            #print("nt : "+str(nt))
                            if(nt["score"] > hiScore):
                                nTs_Hi.append(nt)
                            else: 
                                if(nt["score"] < lowScore and nt["score"] != -999):
                                    nTs_Low.append(nt)
                                else:
                                    if(nt["score"] != -999):
                                        nTs_Bg.append(nt)
                            
                            
                        
                        
                        (nodeTr,linkTr,maxMcff,minMcff,maxSO,minSO) = self.mAdToD3Graph(self.sESTabMcff,self.sESTabsubOpt,self.adjListMcff,self.adjListSO)
                        print("\n\n\ngt\n\n\n")
                        graph_transition = {}
                        graph_transition["nodes"] = nodeTr
                        graph_transition["links"] = linkTr
                        graph_transition["maxMcff"] = maxMcff
                        graph_transition["minMcff"] = minMcff
                        graph_transition["maxSO"] = maxSO
                        graph_transition["minSO"] = minSO
                        
                        tabScoreFiltered = self.filterMinus999(e["scoreTab"])
                        
                        scoresMoy = statistics.mean(tabScoreFiltered)
                        scoreEcart_type =  statistics.stdev(tabScoreFiltered)
                        
                        self.RNA_tab[rna_id_num]["info"] = self.d
                        self.RNA_tab[rna_id_num]["rna_id"] = self.ms_id
                        self.RNA_tab[rna_id_num]["nts"] = nTs
                        self.RNA_tab[rna_id_num]["ScoresMoy"] = scoresMoy
                        self.RNA_tab[rna_id_num]["frequenceNT"] = {}
                        self.RNA_tab[rna_id_num]["frequenceNT"]["A"] = self.calculFreqR("A",seq)
                        self.RNA_tab[rna_id_num]["frequenceNT"]["C"] = self.calculFreqR("C",seq)
                        self.RNA_tab[rna_id_num]["frequenceNT"]["G"] = self.calculFreqR("G",seq)
                        self.RNA_tab[rna_id_num]["frequenceNT"]["U"] = self.calculFreqR("U",seq)
                        self.RNA_tab[rna_id_num]["ncmFreqGlobal"] = self.addD(self.ncmD,{})
                        self.RNA_tab[rna_id_num]["ncmFreqGlobal_mfe"] = self.addD(self.ncmD_mfe,{})
                        self.RNA_tab[rna_id_num]["ncmFreqGlobal_detail"] = self.addD(self.ncmD_detail,{})
                        self.RNA_tab[rna_id_num]["ncmFreqGlobal_mfe_detail"] = self.addD(self.ncmD_mfe_detail,{})
                        self.RNA_tab[rna_id_num]["t1_mcff"] = t1_mcff
                        self.RNA_tab[rna_id_num]["t1_so"] = t1_so
                        self.RNA_tab[rna_id_num]["lowScore"] = lowScore
                        self.RNA_tab[rna_id_num]["hiScore"] = hiScore
                        #self.RNA_tab[rna_id_num]["structure"] = e["structure"]
                        self.RNA_tab[rna_id_num]["stn"] = e["stn"]
                        self.RNA_tab[rna_id_num]["nSoSe_so"] = self.nSoSe_so
                        self.RNA_tab[rna_id_num]["nSoSe_mcff"] = self.nSoSe_mcff
                        self.RNA_tab[rna_id_num]["frP_mcff"] = frP_mcff
                        self.RNA_tab[rna_id_num]["frP_so"] = frP_so
                        self.RNA_tab[rna_id_num]["scoreTab"] = e["scoreTab"]
                        self.RNA_tab[rna_id_num]["erreurTab"] = e["errorTab"]
                        self.RNA_tab[rna_id_num]["seq"] = e["seq"]
                        self.RNA_tab[rna_id_num]["ed"] = self.ed
                        self.RNA_tab[rna_id_num]["freqMfe"] = self.freqMfe
                        self.RNA_tab[rna_id_num]["graph_transition"] = graph_transition
                        self.RNA_tab[rna_id_num]["dotBr2p"] = dotBr2p
                        self.RNA_tab[rna_id_num]["d3ForceLayout2p"] = d3ForceLayout2p
                        self.RNA_tab[rna_id_num]["_2NtMotif"] = _2NtMotif
                        self.RNA_tab[rna_id_num]["scoreEcart_type"] = scoreEcart_type
                        self.RNA_tab[rna_id_num]["TSV_ID"] = e["nid"] 
                        self.RNA_tab[rna_id_num]["reactivityVector"] = self.reactivityVector
                        self.RNA_tab[rna_id_num]["sc_so"] = [sum(x) for x in sc_so]
                        self.RNA_tab[rna_id_num]["sc_mcff"] = [sum(x) for x in sc_mcff]
                        
                        print("stn : "+str(e["stn"]))
                        print("file: "+self.publicFolder+self.rnaJSONFolder+"/"+str(self.ms_id)+".json")
                        f = open(self.publicFolder+self.rnaJSONFolder+"/"+str(self.ms_id)+".json",'w')
                        f.write(json.dumps(self.RNA_tab[rna_id_num]))
                        f.close()
                        
                        for d in dataSetTab:
                            self.newDS(d,self.RNA_tab[rna_id_num])
                        
                        #except:
                          #print("probleme with :"+ str(scoreT)+"\nfile:"+str(name))
            #self.fncm.close()
          else:
            for d in dataSetTab:
              self.newDS(d,json.loads(open(self.publicFolder+self.rnaJSONFolder+"/"+str(self.ms_id)+".json",'r').read()))
          #self.db.close()
      
      def loadCaracteristic(self):
          self.so_ncmCTab = open("listC/so_ncm.txt","r").read().split("\n")
          self.mcff_ncmCTab = open("listC/mcff_ncm.txt","r").read().split("\n")
          self.so_detailCTab = open("listC/so_detail.txt","r").read().split("\n")
          self.mcff_detailCTab = open("listC/mcff_detail.txt","r").read().split("\n")
          
          
          self.listC = {}
          
          self.listC["so_ncm"] = self.so_ncmCTab
          self.listC["mcff_ncm"] = self.mcff_ncmCTab
          self.listC["so_detail"] = self.so_detailCTab
          self.listC["mcff_detail"] = self.mcff_detailCTab
          
          
      def tabToD(self,t):
          d = {}
          for e in t:
            d[e] = 1
          return d
          
          
      
      
      def createNt(self,id_str,voisinG,erreur,score,sorte,position,etat,localNcmD,localNcmD_mfe,localNcmD_detail,localNcmD_mfe_detail,rna_id,inHeliceTemp,freqPairee_mcff,freqPairee_so,seqMotif3,seqMotif5,mfp_mcff,mfp_so,region,voisinPairAllsub_so,voisinPairAllsub_mcff,cpt_so,cpt_mcff,reactivity_pred,ncmTabDG_so_i,ncmTabDG_mcff_i,sc_mcff,sc_so,folder):
          #Voisin droit sera completez apres
          globalRNAinfo = {}
          d_temp = {}
          d_temp["u_id"] = id_str,
          d_temp["voisinGauche"] = voisinG
          d_temp["erreur"] = erreur
          d_temp["score"] = score
          d_temp["sorte"] = sorte
          d_temp["position"] = position
          d_temp["etat"] = etat
          d_temp["localNcmD"] = localNcmD
          d_temp["localNcmD_mcff"] = self.toToWL(localNcmD["mcff"])
          d_temp["localNcmD_so"] = self.toToWL(localNcmD["so"])
          d_temp["localNcmD_mfe"] = localNcmD_mfe
          d_temp["localNcmD_mfe_mcff"] = self.toToWL(localNcmD_mfe["mcff"])
          d_temp["localNcmD_mfe_so"] = self.toToWL(localNcmD_mfe["so"])
          d_temp["localNcmD_detail"] = localNcmD_detail
          d_temp["localNcmD_detail_so"] = self.toToWL(localNcmD_detail["so"])
          d_temp["localNcmD_detail_mcff"] = self.toToWL(localNcmD_detail["mcff"])
          d_temp["localNcmD_mfe_detail"] = localNcmD_mfe_detail
          d_temp["localNcmD_mfe_detail_so"] = self.toToWL(localNcmD_mfe_detail["so"])
          d_temp["localNcmD_mfe_detail_mcff"] = self.toToWL(localNcmD_mfe_detail["mcff"])
          d_temp["cpt_so"] = cpt_so
          d_temp["cpt_mcff"] = cpt_mcff
          d_temp["folder"] = folder
          d_temp["rna_id"] = rna_id
          d_temp["inHelice"] = inHeliceTemp
          d_temp["freqPairee_mcff"] = freqPairee_mcff
          d_temp["freqPairee_so"] = freqPairee_so
          d_temp["seqMotif3"] = seqMotif3
          d_temp["seqMotif5"] = seqMotif5
          d_temp["mfp_mcff"] = mfp_mcff
          d_temp["mfp_so"] = mfp_so
          d_temp["region"] = region
          d_temp["voisinPairAllsub_so"] = voisinPairAllsub_so
          d_temp["voisinPairAllsub_mcff"] = voisinPairAllsub_mcff
          d_temp["reactivity_pred"] = reactivity_pred
          d_temp["ncmTabDG_so"] = ncmTabDG_so_i
          d_temp["ncmTabDG_mcff"] = ncmTabDG_mcff_i
          d_temp["sc_mcff"] = sc_mcff
          d_temp["sc_so"] = sc_so
          
          if(len(sc_mcff) > 0):
            d_temp["sc_sum_mcff"] = sum(sc_mcff)/len(sc_mcff)
          else:
            d_temp["sc_sum_mcff"] = 0 
            
          if(len(sc_so) > 0):
            d_temp["sc_sum_so"] = sum(sc_so)/len(sc_so)
          else:
            d_temp["sc_sum_so"] = 0  
          
          
          return d_temp
          
      def localNcmD_gen(self,mcffNts,soNts,i):
          localNcmD = {}
          localNcmD["mcff"] = mcffNts[i]
          localNcmD["so"] = soNts[i]
          localNcmD["mcff_so"] = extend3(localNcmD["so"],localNcmD["mcff"])
          localNcmD["mcff"] = normalise(localNcmD["mcff"],"sum")#can be max (not implemented)
          localNcmD["so"] = normalise(localNcmD["so"],"sum")
          localNcmD["mcff_so"] = normalise(localNcmD["mcff_so"],"sum")
          return localNcmD
      
      def newDS(self,d,rna):
          tabToWrite = []
          fileName = d["name"]
          hi_threshold = rna["hiScore"]
          low_threshold = rna["lowScore"]
          #print("rna : ",str(rna.keys()))
          #print("len(rna[nts])"+str(len(rna["nts"])))
          #pour tous les nucleotides
          for i in range(0,len(rna["nts"])):
            tabToWrite.append([])
            self.u_id = str(int(self.u_id) + 1)
            score = float(rna["scoreTab"][i])
            #print("d[path] : ",str(d["path"]))
            #print("name : ",str(d["name"]))
            for path in d["path"]:
                ##print("path : ",str(path))
                if(path == 'u_id'):
                  #print("self.u_id : "+str(self.u_id))
                  tabToWrite[-1].append(str(self.u_id))
                else:
                  if(path == 'scoreLabel'):
                    if(score > hi_threshold):
                      tabToWrite[-1].append("Hi")
                    else: 
                      if(score < low_threshold and score != -999):
                        tabToWrite[-1].append("Low")
                      else:
                        if(score != -999):
                          #tabToWrite[-1].append("Bg")
                          tabToWrite.pop()
                          break
                        else:
                          tabToWrite.pop()
                          break
                  else:
                    v = rna
                    #print("path: "+str(path))
                    for pa in path:
                      #print("pa : "+str(pa))
                      #print("type(v) : "+str(type(v)))
                      #if(type(v) == type({})):
                        #print("v : "+str(v.keys()))
                      if(pa == "position"):
                        v = i
                      else:
                        if(pa == "i"):
                          #print("v len : "+str(len(v)))
                          v = v[i]
                          #print("v : "+str(v))
                        else:
                          if(v == None):
                            v = "-"
                            break
                          else:
                            v = v[pa]
                    if(type(v) == type({})):
                      self.addNcm(d["listC"],tabToWrite[-1],v)
                    else:
                      tabToWrite[-1].append(v)
          
          if(self.printCsv):
            f = open(self.publicFolder+self.csvFolder+"/"+self.expName+"_"+fileName+".csv",'a')
            f.write('\n'+'\n'.join([",".join([json.dumps(x).replace(",","||") for x in tab]) for tab in tabToWrite]))
            f.close()  
            
            
            
          return tabToWrite

      def addNcm(self,name,tab,d):
        listCTab = self.listC[name]
        for c in listCTab:
          if(c in d):
            v = d[c]
          else:
            v = 0
          tab.append(v)
      
      
      
      def findMostfreqPartner(self,tabPairTab):
          t1 = []
          if(len(self.seq) > 0):
            for i in range(0,len(self.seq)):
              t1.append([0]*len(self.seq))
            for pT in tabPairTab:
              for p in pT:
                #print("p : "+str(p))
                #print("len : "+str(len(tabPairTab)))
                if(len(p)>1):
                  t1[p[0]][p[1]] += 1
                  t1[p[1]][p[0]] += 1
                else:
                  break
          t2 = []
          for nt in t1:
            pairP = self.compileFreq(nt)
            #print("pairP : "+str(pairP))
            t2.append(pairP)
            
          return (t1,t2) 
              
      def compileFreq(self,tab):
        maxi = 0
        newK = -1
        for i in range(0,len(tab)):
          #print("tab[i] : "+str(tab[i]))
          if(tab[i] > maxi):
            newK = i
            maxi = tab[i]
        if(maxi == 0 ):
          newK = -1
          
        return newK
        
        
        
      def convertDBTo01(self,c):
        if(c == '.'):
          return 0
        else:
          return 1
        
      
      def sumNtWise(self,tabOfStruct):
        t = []
        if(len(tabOfStruct) > 0):
          t = [0]*len(tabOfStruct[0])
          for s in tabOfStruct:
            for i in range(0,len(s)):
              t[i] += self.convertDBTo01(s[i])
          for i in range(0,len(t)):
            t[i] = t[i]/len(tabOfStruct)
        return t
          
      def defVar(self,d,linesParsed):
          self.d = d
          self.testNCM = 4
          self.filePath  = d["filePath"]
          self.expName = d["currentExp"]
          self.folder = d["folder_in"]
          self.inputFileLoader = d["inputFileLoader"]
          self.name = ""+str(self.expName)+"_"+self.folder+"_"+self.inputFileLoader+""
          self.seqProfilTab = linesParsed
          self.RNA_tab = [{}]*len(self.seqProfilTab)
          self.RNA_tab_Str = ''
          
          self.tabSubOpt_allRNA_so = []
          self.tabSubOpt_allRNA_mcff = []
          
          self.indexTabAdded = []
          self.mcffTreshold = d["mcffTreshold"]
          self.soTreshold = d["soTreshold"]
        
          #pour donner une id unique a chaque sous-optimaux (structure secondaire)
          self.indexTabSO = {}
          self.indexTabMcff = {}
          
          #p = pca.PCA( array, fraction=fraction )
      
          self.ms_id = d["uniqueId"]
          self.u_id = d["uniqueId"]
          self.publicFolder = d["publicFolder"]
          self.rnaJSONFolder = d["RNAJSON"]
          self.rnaJSONFolder_in = d["RNAJSON_in"]
          self.csvFolder = d["CSVFolder"]
          self.prediction = d["prediction"]
          self.printCsv = d["printCsv"]
          
          path = self.publicFolder+self.rnaJSONFolder
          if not os.path.exists(path):
                            os.makedirs(path)
                            
          path = self.publicFolder+self.csvFolder
          if not os.path.exists(path):
                            os.makedirs(path)
      
      
      def filterMinus999(self,tab):
        newTab = []
        for x in tab:
          if(x != -999):
            newTab.append(x)
        return newTab
      
      
      def _parse_annotations(self, s):
          d = {}
          s = s.split("|")
          token = ':'
          for item in s:
              if item:
                  pair = item.split(token)
                  if pair[0].strip() in d:
                      d[pair[0].strip()].append(':'.join(pair[1:]))
                  else:
                      d[pair[0].strip()] = [':'.join(pair[1:])]
          return d
          
      def eraseFileAndFillHeader(self,dTab):
          if not os.path.exists(self.publicFolder+self.csvFolder+"/"):
                          os.makedirs(self.publicFolder+self.csvFolder+"/")
          fn = self.publicFolder+self.csvFolder+"/"+self.name+"_Hi"+".csv"
          os.remove(fn) if os.path.exists(fn) else None
          
          fn = self.publicFolder+self.csvFolder+"/"+self.name+"_low"+".csv"
          os.remove(fn) if os.path.exists(fn) else None
          
          fn = self.publicFolder+self.csvFolder+"/"+self.name+"_bg"+".csv"
          os.remove(fn) if os.path.exists(fn) else None
          
          for d in dTab:
            #print("keys : "+str(d.keys()))
            fn = self.publicFolder+self.csvFolder+"/"+self.expName+"_"+d["name"]+".csv"
            os.remove(fn) if os.path.exists(fn) else None
            f = open(fn,'w')
            f.write(",".join(d["header"]) )
            f.close()
      
      def pcaFromCsv(self,fileName):
          
          return p
      
      
      def __str__(self):
          return "self.seqTab : "+",".join(triplet[0] for triplet in self.seqProfilTab)
      
      def findRegion(self,nt_i,l):
          if(nt_i/l < 0.05):
            return "_5p"
          if(nt_i/l > 0.05):
            return "_3p"
          return "millieu"
      


      
          
      def createVoisin(self,voisin):
          d_temp = {}
          d_temp["u_id"] = voisin["u_id"]
          d_temp["etat"] = voisin["etat"]
          d_temp["erreur"] = voisin["erreur"]
          d_temp["score"] = voisin["score"]
          d_temp["sorte"] = voisin["sorte"]
          return d_temp
      
      
      def calculFreqR(self,nt,seq):
          seq = seq.upper()
          freqR = seq.count(nt) / float(len(seq))
          return freqR
      
      def extractMotif(self,n,seq):
          alphabets = ['A', 'C', 'U', 'G']
          keywords = itertools.product(alphabets, repeat = n)
          tab = ["".join(x) for x in keywords ]
          d_temp = {}
          for e in tab:
              d_temp[e] = seq.count(e)
          return d_temp
      
      
      def printinfo(self):
          print ("name : "+self.name)
      
      def jsNodeAll(self,tab,n):
          return [self.jsNode(pair[0],pair[1],1) for pair in tab]+self.phosphoLink(n)
      
      def jsNode(self,source,target,value):
          d = {}
          d["source"] = source
          d["target"] = target
          d["value"] = value
          return d

      def phosphoLink(self,n):
          tab = [{}]*(n-1)
          for i in range(0,n-1):
                  tab[i] = {}
                  tab[i]["source"] = i
                  tab[i]["target"] = i+1
                  tab[i]["value"] = 2
          return tab

      def mcff(self,seq,adjMat):
          
          if (1):#(len(seq) < self.mcffTreshold):
                  self.isThereAnAux = 1
                  stri = 'mcff'+str(" -s " + "'"+seq.upper()+ "' -ft "+str(self.mcffTreshold)+" ")
                  print ("mcffresult : "+stri)
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate()
                  lineSplited = output.strip().decode('ascii').split()
                  print("len lineSplited mcff:"+str(len(lineSplited)))
                  self.sESTabMcff = self.parseMcff(lineSplited)
                  #structureEnergieTab = structureEnergieTab[:80]
                  self.findPairMcff()
                  if(1):#(len(seq) < self.mcffTreshold):
                      self.adjMatMcff = self.adMat2(self.sESTabMcff)
                  self.adjListMcff = self.createAdList(self.indexTabMcff,self.adjMatMcff)

          return self.sESTabMcff


      def parseMcff(self,lineSplited):
          structureEnergyTab = []
          rank = 0
          #self.sESTabMcff = []
          
          for index in range(0,len(lineSplited)):
                  if(lineSplited[index] == "bad"):
                      structureEnergyTab = []
                      break
                  #construction du tableau de sous-optimaaux
                  i=index%3
                  if i == 0:
                      #print ("lineSplited[index] : "+str(lineSplited[index]))
                      struct = lineSplited[index]
                      #if ('(' not in struct):
                              #break
                  elif i== 1 :
                      energy = lineSplited[index]
                  elif i == 2:
                      #print(energy)
                      summary = lineSplited[index]
                      sE = StructEnergy(struct,self.seq,float(energy),summary,rank)
                      #print "strucEnergy McFF : " + sE.displayStructEnergy()
                      structureEnergyTab.append(sE)
                      self.indexTabMcff[struct] = math.floor(index/3)
                      rank = rank+1
          self.nSoSe_mcff = len(structureEnergyTab) / len(self.seq)
          return structureEnergyTab[:self.mcffTreshold]

      def rnaSO(self,seq,adjMat):
          if (1):#len(seq)<self.soTreshold:
                  
                  
                  structureEnergieTab = []
                  #r=10
                  #self.currentrange = (r*self.step,r+self.numStep*self.step)
                  #print "range: ["+str(r*self.step)+','+str(self.smallestTab[r])+']'
                  self.stableRTab = []
                  self.stableCutedTab = []
                  if not os.path.exists(self.filePath+"fa/"):
                          os.makedirs(self.filePath+"fa/")
                  if not os.path.exists(self.filePath+"fa/"+self.folder):
                          os.makedirs(self.filePath+"fa/"+self.folder)
                  fastaP = self.filePath +"fa/"+self.folder+"/" + self.name+".fa"
                  #print ("fastaP = "+fastaP)
                  self.fastaWrite(fastaP,seq)
                  #self.namew = self.name+str(r)
                  stri = "RNAsubopt -e "+str(self.d["so_e_value"])+" -s < \""+fastaP+"\""
                  print (stri)
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate()
                  #print ("out : "+output.strip().decode('ascii'))
                  lineSplited = output.strip().decode('ascii').split()
                  #la sequence
                  #print lineSplited[3]
                  #print ("line de SO : "+str(len(lineSplited)))
                  if(len(lineSplited)>3):
                          self.seq = lineSplited[3]
                          structureEnergieTab = self.parseSO(lineSplited[6:])
                          #print "ParseSo"
                          #structureEnergieTab = structureEnergieTab[:40]
                          #print("self.sESTabsubOpt: "+ "|".join([str(x) for x in self.sESTabsubOpt]))
                          print("structureEnergieTab len : "+str(len(structureEnergieTab)))
                          if(1):#(len(seq) < self.mcffTreshold):
                              self.adjMatSO = self.adMat2(structureEnergieTab)
                          self.createAdjListSO(self.filePath+"Mam")
                          #printPair(self.pairTabSo)

                          #----------------------
                          #GraphMaker_compare(self.seqTab,self.pairTabTab,self.name)#-----------------------

                          #print ("stable region : "+str(self.stableR))
                  return structureEnergieTab
        

      def parseSO(self,lineSplited):
          finalTab = []
          rnaSuboptTab = []
          for i in range(0,len(lineSplited)):
                      #print("i : "+str(i)+" line : "+str(lineSplited[i]))
                      if ( i%2 == 1):
                          #print("1 : " + lineSplited[i])
                          e = lineSplited[i]
                          rnaSuboptTab.append(st+" "+str(e))
                      else:
                          #print("0 : " + lineSplited[i])
                          st = lineSplited[i]
          print ("nombre de SO : "+str(len(rnaSuboptTab)))
          if len(rnaSuboptTab)>0:
                      #print "mfe_SO : "+rnaSuboptTab[0]
                      #print "#SO : "+str(len(rnaSuboptTab))
                      
                      for i in range(0,len(rnaSuboptTab)):
                          #print ("rnaSuboptTab["+str(i)+"] : "+rnaSuboptTab[i])
                          rsoSplited = rnaSuboptTab[i].split()
                          if (float(rsoSplited[1]) == 0):
                                  #print("breakOn : "+rsoSplited[0])
                                  break
                          sE = StructEnergy(rsoSplited[0],self.seq,float(rsoSplited[1]),"-",i)
                          #print("####rsoSplited[0] : "+rsoSplited[0])
                          self.indexTabSO[rsoSplited[0]] = i
                          #print ("strucEnergy SO : " + sE.displayStructEnergy())
                          finalTab.append(sE)
          self.nSoSe_so = len(finalTab) / len(self.seq)
          return finalTab[:self.soTreshold]
        

      def rnaFold(self,seq,adjMat):
          if (1):#len(seq)<self.soTreshold:
                  
                  
                  structureEnergieTab = []
                  #r=10
                  #self.currentrange = (r*self.step,r+self.numStep*self.step)
                  #print "range: ["+str(r*self.step)+','+str(self.smallestTab[r])+']'
                  self.stableRTab = []
                  self.stableCutedTab = []
                  if not os.path.exists(self.filePath+"fa/"):
                          os.makedirs(self.filePath+"fa/")
                  if not os.path.exists(self.filePath+"fa/"+self.folder):
                          os.makedirs(self.filePath+"fa/"+self.folder)
                  fastaP = self.filePath +"fa/"+self.folder+"/" + self.name+".fa"
                  #print ("fastaP = "+fastaP)
                  self.fastaWrite(fastaP,seq)
                  #self.namew = self.name+str(r)
                  stri = "RNAfold -p -d2 --noLP −−noPS < \""+fastaP+"\""
                  #print (stri)
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate()
                  #print ("out : "+output.strip().decode('ascii'))
                  lineSplited = output.strip().decode('ascii').split()
                  #la sequence
                  #print lineSplited[3]
                  print ("line de Fold : "+str(len(lineSplited)))
                  if(len(lineSplited)>0):
                          seqFold = lineSplited[1]
                          (ed,freqMfe) = self.parseFold(lineSplited[2:])
                          print("Ensemble diversity : "+str(ed)+" | frequence Mfe : "+str(freqMfe))
                          
                          #print("structureEnergieTab len : "+str(len(structureEnergieTab)))
                          #if(1):#(len(seq) < self.mcffTreshold):
                              #self.adjMatSO = self.adMat2(structureEnergieTab)
                          #self.createAdjListSO(self.filePath+"Mam")
                          #printPair(self.pairTabSo)

                          #----------------------
                          #GraphMaker_compare(self.seqTab,self.pairTabTab,self.name)#-----------------------

                          #print ("stable region : "+str(self.stableR))
                  return (ed,freqMfe)
        

      def parseFold(self,lineSplited):
          finalTab = []
          
          freqMfe = -999
          ed = -999
          
          nextIsfreqMfe = False
          nextIsde = False
          for i in range(0,len(lineSplited)):
                      #print("i : "+str(i)+" line : "+str(lineSplited[i]))
                      if(nextIsfreqMfe):
                        freqMfe = float(lineSplited[i][:-1])
                        nextIsfreqMfe = False
                      if(lineSplited[i] == 'ensemble' and freqMfe == -999):
                        nextIsfreqMfe = True
                      
                      if(nextIsde):
                        ed = float(lineSplited[i])
                        nextIsde = False
                      if(lineSplited[i] == 'diversity'):
                        nextIsde = True
                      
          return (ed,freqMfe)

      def parseFasta(self,fastaStr):
          seqTab = []
          s = ""
          n = ""
          fTab = fastaStr.split("\n")
          for line in fTab:
                      if(line.startswith( '>' )):
                          if(s != ""):
                                  seqTab.append({'name' : n,'seq':s}) 
                                  n = line[1:]
                                  s = ""
                          else:
                                  n = line[1:]
                      else:
                          s += line.rstrip()
          seqTab.append({'name' : n,'seq': s}) 
          return seqTab

      def fastaWrite(self, outFile,seq):
        if not os.path.exists(outFile+"fa/"):
            try:  
                    with open(outFile,"w") as newFasta:
                            header = ">" + self.name + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                            newFasta.write(header + seq)
            except IOError:
                    print ("Failed to open " + outFile)
                    exit(1)

        #print ">> Done!"
      
      def linkForce(self,rna1,rna2):
          res = -1
          alignement = pairwise2.align.globalxx(rna1,rna2)
          longuest = max(len(rna1),len(rna2))
          if (len(alignement) < 1):
                  print("Alignement error")
                  res = -1
          else:
                  res = alignement[0][2]/longuest
          return res
      
      
      
      def createRNALink(self):
          returnTab = []
          for i_rna in range(0,len(self.seqProfilTab)):
                  print("Comparing seq : " +self.seqProfilTab[i_rna]["seq"])
                  for j_rna in range(i_rna+1,len(self.seqProfilTab)):
                          f = self.linkForce(self.seqProfilTab[i_rna]["seq"],self.seqProfilTab[j_rna]["seq"])
                          if(f>0.9):
                                  returnTab.append("{\"source\":"+str(i_rna)+",\"target\":"+str(j_rna)+",\"value\":"+str(f)+"}")
          return returnTab
      
      
      def findPairSo(self):
          print ("startFindPair")
          if(len(self.sESTabsubOpt)>0):
                  #print self.sESTabsubOpt[0].struct
                  stList = list(self.sESTabsubOpt[0].struct)
                  #Pour tous les nucleotides (, trouve sa paire
                  self.pairTabSo = []
                  for n in range(0,len(stList)):
                          if stList[n] == '(':
                                  i = self.findPartnerFromStruct(n,self.sESTabsubOpt[0].struct)
                                  if(i == -1):
                                          print ("pairNOTfoundSO")
                                          break
                                  self.pairTabSo.append((n,i))
                  print ("pairfoundSO")

      def findPairMcff(self):
          if(len(self.sESTabMcff)>0):
                  #print self.SEnergyTabMcff[0].struct
                  stList = list(self.sESTabMcff[0].struct)
                  #Pour tous les nucleotides (, trouve sa paire
                  for n in range(0,len(stList)):
                          if stList[n] == '(':
                                  i = self.findPartnerFromStruct(n,self.sESTabMcff[0].struct)
                                  self.pairTabMcff.append((n,i)) 
                  print ("pairfound")


      def findPairGen(self,secStruct):
          openDict = {}
          openDict['('] = 1 
          openDict['['] = 1 
          openDict['{'] = 1 
          openDict['<'] = 1 
          openDict['A'] = 1 
          openDict['B'] = 1 
          tempTab=[]
          #print self.SEnergyTabMcff[0].struct
          stList = list(secStruct)
          #Pour tous les nucleotides (, trouve sa paire
          for n in range(0,len(stList)):
                  c = stList[n]
                  if c in openDict :
                          i = self.findPartnerGen(n,secStruct)
                          if(i == -1):
                                  print ("findPAirNOMA : " + str(n))
                          tempTab.append((n,i))
          return tempTab
          
        
      def voisinDroit(self,i,l,d):
          if(i < l):
              #print("vg i : "+str(i))
              for j in range(i+1,l):
                  #print("j : "+str(j))
                  if(j in d):
                      
                      return j
              return -1
          else:
              return -2

      def voisinGauche(self,i,d):
          if(i > 0):
              #print("vd i : "+str(i))
              for j in range(i-1,-1,-1):
                  #print("j : "+str(j))
                  if(j in d):
                      return j
              return -1
          else:
              return -2
  
        
      def pbToNcmAll(self,l,detailed):
              
              l_so = len(self.pairTabAllSubOpt_so)
              l_mcff = len(self.pairTabAllSubOpt_mcff)
              finalD = {}
              finalD['mcff'] = [{}]*l_mcff
              finalD['mcff_merged'] = {}
              finalD['so'] = [{}]*l_so
              finalD['so_merged'] = {}
              colorPredTab_so = [[] for x in range(0,l)]
              soNts = [{}]*l
              self.tabVoisinAllsub_so = []
              ncm_so = []
              #pour tous les sous-optimaux de RNASubOpt
              sCTab_by_subOpt_so = []
              for j in range(0,len(self.pairTabAllSubOpt_so)): #pour tous les sous-optimaux (so)
                  so = self.pairTabAllSubOpt_so[j]
                  #f.write(";".join("("+str(x[0])+"|"+str(x[1])+")" for x in so)+"\n")
                  (tabNtsD,tabVoisin1,ncmPredictionNtTab_so,ncm_temp,scoreCoherenceTab_so) = self.pbToNcm(so,l,"so",detailed)
                  sCTab_by_subOpt_so.append(scoreCoherenceTab_so)
                  self.tabVoisinAllsub_so.append(tabVoisin1)
                  ncm_so.append(ncm_temp)
                  #print("ncm_temp : "+str(ncm_temp))
                  #print("tabNtsD : "+";".join(str(x) for x in tabNtsD))
                  #print("d : "+str(tabNtsD))
                  for i in range(0,len(tabNtsD)):
                      #print("colorPredTab_so : " + str(len(colorPredTab_so)))
                      #print("i : " + str(i))
                      if(len(ncmPredictionNtTab_so) > 0):
                        colorPredTab_so[i].append(ncmPredictionNtTab_so[i])
                      d = {}
                      d.clear()
                      d = tabNtsD[i]
                      soNts[i] = self.addD(soNts[i],d)
                      finalD['so'][j] = self.addD(finalD['so'][j],d)
                      finalD['so_merged'] = self.addD2(finalD['so_merged'],d)
                      #print("finalD['so_merged']:"+str(len(finalD['so_merged'])))
                      #finalD['so_merged'] = extend_so(finalD['so_merged'],self.ncmCaracteristic_d)
              #print("len(finalD['so_merged'].keys()) : "+str(len(finalD['so_merged'].keys())))
              self.tabVoisinAllsub_mcff = []
              colorPredTab_mcff = [[] for x in range(0,l)]
              mcffNts = [{}]*l
              ncm_mcff = []
              sCTab_by_subOpt_mcff = []
              for j in range(0,len(self.pairTabAllSubOpt_mcff)):#pour tous les sous-optimaux (mcff)
                  mcff = self.pairTabAllSubOpt_mcff[j]
                  (tabNtsD,tabVoisin2,ncmPredictionNtTab_mcff,ncm_temp,scoreCoherenceTab_mcff)  = self.pbToNcm(mcff,l,"mcff",detailed)
                  sCTab_by_subOpt_mcff.append(scoreCoherenceTab_mcff)
                  self.tabVoisinAllsub_mcff.append(tabVoisin2)
                  ncm_mcff.append(ncm_temp)
                  #print("ncm_temp : "+str(ncm_temp))
                  #print("tabNtsD : "+";".join(str(x) for x in tabNtsD))
                  for i in range(0,len(tabNtsD)):
                      if(len(ncmPredictionNtTab_mcff) > 0):
                        colorPredTab_mcff[i].append(ncmPredictionNtTab_mcff[i])
                      d = {}
                      d.clear()
                      d = tabNtsD[i]
                      mcffNts[i] = self.addD(mcffNts[i],d)
                      finalD['mcff'][j] = self.addD(finalD['mcff'][j],d)
                      finalD['mcff_merged'] = self.addD2(finalD['mcff_merged'],d)
                      #print("finalD['mcff_merged']"+str(finalD['mcff_merged']))
                      #finalD['mcff_merged'] = extend_mcff(finalD['mcff_merged'],self.ncmCaracteristic_d)
              #print("len(finalD['mcff_merged'].keys()) : "+str(len(finalD['mcff_merged'].keys())))
              #print("finalD['mcffNts'] : "+ str(finalD['mcffNts']))
              #print("finalD['soNts'] : "+ str(finalD['soNts']))
              #print(finalD['mcffNts'])
              #print("ncm_so : " +str(ncm_so))
              #print("ncm_mcff : " +str(ncm_mcff))
              ncm_so_t = self.transpose(ncm_so)
              ncm_mcff_t = self.transpose(ncm_mcff)
              sc_so = sCTab_by_subOpt_so
              sc_mcff = sCTab_by_subOpt_mcff
              return (finalD,mcffNts,soNts,colorPredTab_so,colorPredTab_mcff,ncm_so_t,ncm_mcff_t,sc_so,sc_mcff)
      
      def transposeSc(self,tabt):
        tab = [ [] for x in tabt[0]]
        #print("len(tabt[0]) : "+str(len(tabt[0])))
        
        for t in tabt:
          for i in range(0,len(t)):
            #print("i : "+str(i))
            #print(str(t))
              
            tab[i].append(t[i])
            
        return tab
      
      def transpose(self,tabt):
        tab = [ [] for x in tabt[0]]
        #print("len(tabt[0]) : "+str(len(tabt[0])))
        
        for t in tabt:
          for i in range(0,len(t)):
            #print("i : "+str(i))
            #print(str(t))
              
            tab[i].append({"droit":t[i]["droit"],"gauche":t[i]["gauche"]})
            
        return tab

      def pbToNcmAll_mfe(self,l,detailed):
              l_so = len(self.pairTabAllSubOpt_so)
              l_mcff = len(self.pairTabAllSubOpt_mcff)
              finalD = {}
              finalD['mcff'] = [{}]*l_mcff
              finalD['mcff_merged'] = {}
              finalD['so'] = [{}]*l_so
              finalD['so_merged'] = {}
              soNts = [{}]*l
              #pour la mfe de RNASubOpt
              for j in range(0,1): #pour tous les sous-optimaux (so)
                  so = self.pairTabAllSubOpt_so[j]
                  #f.write(";".join("("+str(x[0])+"|"+str(x[1])+")" for x in so)+"\n")
                  (tabNtsD,tabVoisin,ncmPredictionNtTab_so,ncm,scoreCoherenceTab_so)  = self.pbToNcm(so,l,"so",detailed)
                  #print("tabNtsD : "+";".join(str(x) for x in tabNtsD))
                  #print("d : "+str(tabNtsD))
                  for i in range(0,len(tabNtsD)):
                      d = {}
                      d.clear()
                      d = tabNtsD[i]
                      soNts[i] = self.addD(soNts[i],d)
                      finalD['so'][j] = self.addD(finalD['so'][j],d)
                      finalD['so_merged'] = self.addD2(finalD['so_merged'],d)
                      #print("finalD['so_merged']:"+str(len(finalD['so_merged'])))
                      #finalD['so_merged'] = extend_so(finalD['so_merged'],self.ncmCaracteristic_d)
              print("len(finalD['so_merged'].keys()) : "+str(len(finalD['so_merged'].keys())))
              mcffNts = l*[{}]
              for j in range(0,1):#pour la mfe (mcff)
                  mcff = self.pairTabAllSubOpt_mcff[j]
                  (tabNtsD,tabVoisin,ncmPredictionNtTab_mcff,ncm,scoreCoherenceTab_mcff)  = self.pbToNcm(mcff,l,"mcff",detailed)
                  #print("tabNtsD : "+";".join(str(x) for x in tabNtsD))
                  for i in range(0,len(tabNtsD)):
                      d = {}
                      d.clear()
                      d = tabNtsD[i]
                      mcffNts[i] = self.addD(mcffNts[i],d)
                      finalD['mcff'][j] = self.addD(finalD['mcff'][j],d)
                      finalD['mcff_merged'] = self.addD2(finalD['mcff_merged'],d)
                      #print("finalD['mcff_merged']"+str(finalD['mcff_merged']))
                      #finalD['mcff_merged'] = extend_mcff(finalD['mcff_merged'],self.ncmCaracteristic_d)
              print("len(finalD['mcff_merged'].keys()) : "+str(len(finalD['mcff_merged'].keys())))
              #print("finalD['mcffNts'] : "+ str(finalD['mcffNts']))
              #print("finalD['soNts'] : "+ str(finalD['soNts']))
              #print(finalD['mcffNts'])
              return (finalD,mcffNts,soNts)             
      
      
      def seqNcm(self,seq1,seq2,detailed):
        if(detailed):
          if(len(seq1) > 6):
            seq1 = "L"
          if(len(seq2) > 6):
            seq2 = "L"
          return "-"+seq1+"-"+seq2
        else:
          return ""
      
      def posNcm(self,pos,detailed):
        if(detailed):
          return "_pos_"+str(pos)
        else:
          return ""
      
      
      def openMongoClient(self):
        d = {}
        client = MongoClient()
        client = MongoClient('10.0.0.161', 27027)
        db = client.rdv
        return db
      
      def getPrediction2(self,low,bg,hi):
        l_low = low
        l_bg = bg
        l_hi = hi
        if((l_low + l_hi + l_bg) > self.d['nbNt_Ncm'] and l_bg / float(l_low + l_hi + l_bg) > self.d["predCutOff"] ):
          return "Bg"
        if((l_low + l_hi + l_bg) > self.d['nbNt_Ncm'] and l_hi / float(l_low + l_hi + l_bg) > self.d["predCutOff"] ):
          return "Hi"
        if((l_low + l_hi + l_bg) > self.d['nbNt_Ncm'] and l_low / float(l_low + l_hi + l_bg) > self.d["predCutOff"] ):
          return "Low"
        
        return "NED"
      
      def getPrediction(self,low,bg,hi):
        l_low = low
        l_bg = bg
        l_hi = hi
        poids = abs(l_hi - l_low)/(l_hi + l_low)
        if((l_low + l_hi ) > self.d['nbNt_Ncm'] and l_hi / float(l_low + l_hi ) > self.d["predCutOff"] ):
          return ("Hi",poids)
        if((l_low + l_hi ) > self.d['nbNt_Ncm'] and l_low / float(l_low + l_hi ) > self.d["predCutOff"] ):
          return ("Low",poids)
        
        return ("NED",0)
      
      def colorTranlation(self,pred):
        d = {"Hi":"blue","Low":"red","NED":"white","Bg":"purple"}
        return d[pred]
      
      def getColor2Ncm(self,p1,p2):
        if(p1 ==  p2):
          return self.colorTranlation(p1)
        if(p1 == "Bg" or p1 == "NED"):
          return self.colorTranlation(p2)
        if(p2 == "Bg" or p2 == "NED"):
          return self.colorTranlation(p1)
        # p1 and p2 are Hi and Low or Low and Hi ...
        return "indigo"
      
      
      def getColorCoherence(self,ncm,i,soft):
        actualLevel = self.reactivityVector[i]
        scoh = 0
        color = "black"
        if(actualLevel == "Bg"):
          color = "grey"
        else:
          if(actualLevel == "NA"):
            color = "black"
          else:
            if(len(ncm) == 2):
              #print("ncm1 : "+ ncm[0])
              cursor_1 = self.db.ncmRdat2_10.find({"ncm":ncm[0],"soft":soft})
              ncm1 = False
              for d1 in cursor_1:
                ncm1 = True
                low1 = d1["low"]
                #print("low1 : "+ str(low1))
                hi1 = d1["hi"]
                #print("hi1 : "+ str(hi1))
                bg1 = d1["bg"]
                #print("bg1 : "+ str(bg1)+"\n")
              
              #print("ncm2 : "+ ncm[1])
              cursor_2 = self.db.ncmRdat2_10.find({"ncm":ncm[1],"soft":soft})
              ncm2 = False
              for d2 in cursor_2:
                ncm2 = True
                low2 = d2["low"]
                #print("low2 : "+ str(low2))
                hi2 = d2["hi"]
                #print("hi2 : "+ str(hi2))
                bg2 = d2["bg"]
                #print("bg2 : "+ str(bg2)+"\n")
                
                
              if(ncm1 and ncm2):
                (pred1,poids) = self.getPrediction(low1,bg1,hi1)
                (pred2,poids) = self.getPrediction(low2,bg2,hi2)
                
                color = self.getColor2Ncm(pred1,pred2)
              else :
                color = "green"
            else:
              #print("ncm : "+ ncm[0])
              cursor = self.db.ncmRdat2_10.find({"ncm":ncm[0],"soft":soft})
              isNcmpresent = False
              for d in cursor:
                isNcmpresent = True
                low = d["low"]
                #print("low : "+ str(low))
                hi = d["hi"]
                #print("hi : "+ str(hi))
                bg = d["bg"]
                #print("bg : "+ str(bg)+"\n")
              if(isNcmpresent):
                (pred,poids) = self.getPrediction(low,bg,hi)
                color = self.getColor2Ncm(pred)
              else :
                color = "green"
          
        if(color == "red" and actualLevel == "Low" or color == "blue" and actualLevel == "Hi"):
          scoh = 1 * poids
        if(color == "blue" and actualLevel == "Low" or color == "red" and actualLevel == "Hi"):
          scoh = -1 * poids
          
        return (color,scoh)
      
      
      
      def pbToNcm(self,paireBaseTab,l_seq,soft,detailed):
              
              #print(",".join("("+str(x[0])+","+str(x[1])+")" for x in paireBaseTab))
              self.db = self.openMongoClient()
              d = {}
              ncmD = []
              
              #print("paireBaseTab: "+";".join("("+str(x[0])+"|"+str(x[1])+")" for x in paireBaseTab))
              for i in range(0,len(paireBaseTab)):
                  #print('value :'+str(paireBaseTab[i]))
                  t = paireBaseTab[i][0]
                  s = paireBaseTab[i][1]
                  d[s] = t
                  d[t] = s
                  #print('s :' + str(s) + ' t = '+str(t))
                  
              tabOfVoisin = l_seq*[None]    
              tabOfNts = l_seq*[None]
              '''
              for key, value in d.items() :
                  print ("pb : "+str(key)+","+str( value))
              '''
              s = ""
              for i in range(0,l_seq):
                  #print("s : "+s)
                  ncmD.append({})
                  ncmD[i].clear()
                  ncmD[i]['droit'] = []
                  ncmD[i]['gauche'] = []
                  tuplet = []
                  tabOfVoisin[i] = {}
                  
                  
                  vd = self.voisinDroit(i,l_seq,d)
                  tabOfVoisin[i]['vd'] = vd
                  if(vd in d):
                    tabOfVoisin[i]['p_vd'] = d[vd]
                  else:
                    tabOfVoisin[i]['p_vd'] = "-"
                  
                  
                  vg = self.voisinGauche(i,d)
                  tabOfVoisin[i]['vg'] = vg
                  if(vg in d):
                    tabOfVoisin[i]['p_vg'] = d[vg]
                  else:
                    tabOfVoisin[i]['p_vg'] = "-"
                    
                  tabOfVoisin[i]['i'] = i  
                  if(vd >= 0 and vg >= 0 and vd == int(d[vg])):
                      if(i not in d):
                        tabOfVoisin[i]['p_i'] = -1
                      else:
                        tabOfVoisin[i]['p_i'] = d[i]
                      #print("--res :"+str(len(range(vd,int(d[vd])+1))))
                      s = str(str(vd-d[vd]+1)+self.seqNcm("-",self.seq[vg:vd+1],detailed)+self.posNcm(str(i-vg),True))
                      ncmD[i]['droit'].append(s)
                      ncmD[i]['gauche'].append("-")
                      tuplet.append(s)
                      tuplet.append(s)#pour balancer
                      #loop
                  else:
                      if(i in d):
                          tabOfVoisin[i]['p_i'] = d[i]
                          if(vd >= 0 ):
                              if(i == d[vd]):
                                  loop_v = vd-i+1
                                  if(loop_v>10):
                                      ncmD[i]['droit'].append("L")
                                      tuplet.append("L")
                                  else:
                                      s = str(str(loop_v)+self.seqNcm("-",self.seq[i-1:vd],detailed)+self.posNcm(str(0),True))
                                      ncmD[i]['droit'].append(s)
                                      tuplet.append(s)
                                  #loop
                              else:
                                  ntBetweVd1 = vd-i+1
                                  ntBetweVd2 = d[i]-d[vd]+1
                                  #print("vg1 : "+ str(vd)+" vg2 : "+ str(d[vd])+" i1 : "+ str(i)+" i2 : "+ str(d[i]))
                                  #print("--res :"+str(ntBetweVd1)+','+str(ntBetweVd2))
                                  if(abs(ntBetweVd1) > 6 ):
                                      ntBetweVd1Str = "L"
                                      seqBt1a = "L"
                                      seqBt1b = "L"
                                  else:
                                      ntBetweVd1Str = str(ntBetweVd1)
                                      seqBt1a = ntBetweVd1Str
                                      seqBt1b = self.seq[i:vd+1]
                                      
                                  if(abs(ntBetweVd2) > 6 ):
                                      ntBetweVd2Str = "L"
                                      seqBt2a = "L"
                                      seqBt2b = "L"
                                  else:
                                      ntBetweVd2Str = str(ntBetweVd2) 
                                      seqBt2a = ntBetweVd2Str
                                      seqBt2b = self.seq[d[vd]:d[i]+1]
                                      
                                  seqBt = seqBt1a+'_'+seqBt2a+self.seqNcm(seqBt1b,seqBt2b,detailed)+self.posNcm(str(0),True)
                                  ncmD[i]['droit'].append(seqBt)
                                  tuplet.append(seqBt)
                              
                          else:
                              if(vd == -1):
                                  ncmD[i]['droit'].append("noVD_Paired")
                                  tuplet.append("noVD_Paired")
                              if(vd == -2):
                                  ncmD[i]['droit'].append("Dernier_Paired")
                                  tuplet.append("Dernier_Paired")
                          #voisin droit   
                          if(vg >= 0 ):
                              if(i == d[vg]):
                                loop_v = (i-vg)+1
                                if(loop_v>10):
                                      loop_v = "L"
                                s = str(str(loop_v)+self.seqNcm("-",self.seq[vg:i+1],detailed)+self.posNcm(str(i-vg),True))
                                ncmD[i]['gauche'].append(s)
                                tuplet.append(s)
                                #loop
                              else:
                                  ntBetweVg1 = i-vg+1
                                  ntBetweVg2 = d[vg]-d[i]+1
                                  #print("vg1 : "+ str(vg)+" vg2 : "+ str(d[vg])+" i1 : "+ str(i)+" i2 : "+ str(d[i]))
                                  #print("--res :"+str(ntBetweVg1)+','+str(ntBetweVg2))
                                  if(ntBetweVg1 > 6 or -1*ntBetweVg1 > 6):
                                      ntBetweVg1Str = "L"
                                  else:
                                      ntBetweVg1Str = str(ntBetweVg1)
                                      
                                  if(ntBetweVg2 > 6 or -1*ntBetweVg2 > 6):
                                      ntBetweVg2Str = "L"
                                  else:
                                      ntBetweVg2Str = str(ntBetweVg2) 
                                      
                                  s = str(ntBetweVg1Str+'_'+ntBetweVg2Str + self.seqNcm(self.seq[vg:i+1],self.seq[d[i]:d[vg]+1],detailed)+self.posNcm(str(i-vg),True))
                                  ncmD[i]['gauche'].append(s)
                                  tuplet.append(s)
                          else:
                              
                              if(vg == -1):
                                  ncmD[i]['gauche'].append("noVG_Paired")
                                  tuplet.append("noVG_Paired")
                              
                              if(vg == -2):
                                  ncmD[i]['gauche'].append("Premier_Paired")
                                  tuplet.append("Premier_Paired")
                          #voisin gauche
                          #deux ncms
                      else:
                          tabOfVoisin[i]['p_i'] = -1
                          if(vd < 0 or vg < 0):
                              if(vd == -1):
                                  ncmD[i]['droit'].append("noVD_notPaired")
                                  ncmD[i]['gauche'].append("-")
                                  tuplet.append("noVD_notPaired")
                                  tuplet.append("noVD_notPaired")
                              if(vg == -1):
                                  ncmD[i]['gauche'].append("noVG_notPaired")
                                  ncmD[i]['droit'].append("-")
                                  tuplet.append("noVG_notPaired")
                                  tuplet.append("noVG_notPaired")
                              if(vd == -2):
                                  ncmD[i]['droit'].append("Dernier_notPaired")
                                  ncmD[i]['gauche'].append("-")
                                  tuplet.append("Dernier_notPaired")
                                  tuplet.append("Dernier_notPaired")
                              if(vg == -2):
                                  ncmD[i]['gauche'].append("Premier_notPaired")
                                  ncmD[i]['droit'].append("-")
                                  tuplet.append("Premier_notPaired")
                                  tuplet.append("Premier_notPaired")
                          else:
                              ntBetweVdEtVg1 = vd-vg+1
                              ntBetweVdEtVg2 = d[vg]-d[vd]+1
                              #print("vg1 : "+ str(vg)+" vg2 : "+ str(d[vg])+" vd1 : "+ str(vd)+" vd2 : "+ str(d[vd]))
                              #print("--res :"+str(ntBetweVdEtVg1)+','+str(ntBetweVdEtVg2))
                              if(ntBetweVdEtVg1 > 6 or -1*ntBetweVdEtVg1 > 6):
                                  ntBetweVdEtVg1Str = "L"
                              else:
                                  ntBetweVdEtVg1Str = str(ntBetweVdEtVg1)
                                  
                              if(ntBetweVdEtVg2 > 6 or -1*ntBetweVdEtVg2 > 6):
                                  ntBetweVdEtVg2Str = "L"
                              else:
                                  ntBetweVdEtVg2Str = str(ntBetweVdEtVg2)    
                              
                              s = str(ntBetweVdEtVg1Str+'_'+ntBetweVdEtVg2Str + self.seqNcm(self.seq[vg:vd+1],self.seq[d[vd]:d[vg]+1],detailed)+self.posNcm(str(i-vg),True))
                              ncmD[i]['droit'].append(s)
                              ncmD[i]['gauche'].append("-")
                              
                              tuplet.append(s)
                              #pour balancer
                              tuplet.append(s)
                              s = ""
                  '''
                  if(soft == "mcff"):
                    tabOfVoisin[i]["i"] = 1
                    tabOfVoisin[i]['p_i'] = 79
                    tabOfVoisin[i]["vg"] = 0
                    tabOfVoisin[i]['p_vg'] = 80
                    tabOfVoisin[i]["vd"] = 2
                    tabOfVoisin[i]['p_vd'] = 78
                  
                  if(i == self.testNCM):
                      print("tuplet : "+str(tuplet)) 
                  
                  if(len(tuplet) != 2):
                      print("tuplet : "+str(tuplet)) 
                  #print("tuplet : "+ str(tuplet)+"\n")
                  '''
                  tabOfNts[i] = tuplet
                  #print("ncmD["+str(i)+"] : "+str(ncmD[i])  )
                  
                  
                  
              ntNcm = [None] * len(tabOfNts)
              ncmPredictionNtTab = []
              scoreTab = []
              for j in range(0,len(tabOfNts)):
                  tup = tabOfNts[j]
                  if(detailed):
                    if(self.prediction):
                      (color,scoh) = self.getColorCoherence(tup,j,soft)
                      #if(win == 'red' or win == 'blue'):
                        #print("color : "+win)
                      ncmPredictionNtTab.append(color)
                      scoreTab.append(scoh)
                    else:
                      ncmPredictionNtTab.append("black")
                      scoreTab.append(1)
                  ntNcm[j] = {}
                  for ncm in tup:
                      #print("ncm : "+ncm)
                      #print("ncm : "+ncm+"_"+soft)
                      if(ncm+"_"+soft in self.cpjp):
                          ncm = "cpjp"
                      
                      #if(ncm+"_"+soft in self.ncmCaracteristic_d):
                      if(ncm in ntNcm[j]):
                          ntNcm[j][ncm] += 1
                      else:
                          ntNcm[j][ncm] = 1
                      #else:
                          #self.fncm.write(str(ncm)+"\n")
                      
                  #ntNcm[j]= extend3(ntNcm[j],self.ncmCaracteristic_d)
              #print("ntNcm.items() : "+";".join(str(x) for x in ntNcm))
              ''' 
              for key, value in ntNcm.items() :
                  print (key, value)
              '''
              #print("ntNcm:"+str(ntNcm))
              #print("ncmD[79] : "+str(ncmD[79]))
              self.db.close
              return (ntNcm,tabOfVoisin,ncmPredictionNtTab,ncmD,scoreTab)
      
      def addD2(self,d1,d2):
          dFin = {}
          for (key,value) in d1.items():
              if key in d2 :
                  dFin[key] = value + d2[key]
              else:
                  dFin[key] = value
          for (key,value) in d2.items():
              if not (key in d1) :
                  dFin[key] = value
          return dFin
          
      def addD(self,d1,d2):
          dFin = {}
          for (key,value) in d1.items():
              if key in d2 :
                  dFin[key] = value + d2[key]
              else:
                  dFin[key] = value
          for (key,value) in d2.items():
              if not (key in d1) :
                  dFin[key] = value
          return dFin
          
          
      def countStableRegion(self,pairNeeded):
              acc = 0
              stablePart = 0
              tempTab = []
              sRFound = False
              for i in range(0,len(self.pairTabSo)-1):
                      if(self.pairTabSo[i][0]+1 == self.pairTabSo[i+1][0] and self.pairTabSo[i][1]-1 == self.pairTabSo[i+1][1]):
                              tempTab.append((self.pairTabSo[i][0],self.pairTabSo[i][1]))
                              acc+=1
                              if (acc == pairNeeded-1):
                                  sRFound = True
                      else:
                              acc = 0
                              if(sRFound):
                                  self.stableRTab.append(tempTab)
                                  stablePart +=1
                              sRFound = False
                              tempTab = []
              if (sRFound):
                      self.stableRTab.append(tempTab)
                      stablePart +=1
              return stablePart

      def appendStableRegion(self):
          if(len(self.sESTabsubOpt)>0):
                  #print "structMFE : "+self.sESTabsubOpt[0].struct
                  for i in range(0,len(self.stableRTab)):
                          cutted = self.cutRNA(i,self.sESTabsubOpt[0].struct)
                          self.seqTab.append(cutted[0])
                          self.pairTabTab.append(cutted[1])
                          reelIndex = self.stableRTab[i][0][0] + self.r*self.step
                          self.stableCutedTab.append(RNA2D_L(self.name+"_Cut_"+str(self.stableRTab[i][0][0]),cutted[0],reelIndex,self.filePath))
                          #break

      def cutRNA(self,i,secStruct):
          firstNt = self.stableRTab[i][0][0]
          lastNt = self.stableRTab[i][0][1]
          lseq = list(self.seq)
          lsecStruct = list(secStruct)
          if(len(self.stableRTab) < i+2 or self.isTerminal(i)):
                  secStructCutted = "".join(lsecStruct[firstNt:lastNt+1])
                  #print "secStructCutted : "+secStructCutted
                  #print "seqCutted : "+"".join(lseq[firstNt:lastNt+1])
                  return ("".join(lseq[firstNt:lastNt+1]),findPairGen(secStructCutted))
          else:
                  firs2tNt = self.stableRTab[i+1][-1][0]
                  last2Nt = self.stableRTab[i+1][-1][1]
                  secStructCutted = "".join(lsecStruct[firstNt:firs2tNt+1])+"...."+"".join(lsecStruct[last2Nt:lastNt+1])
                  seqCuted = "".join(lseq[firstNt:firs2tNt+1]) +"GAAA"+"".join(lseq[last2Nt:lastNt+1])
                  #print "secStructCutted : "+secStructCutted
                  #print "seqCutted : "+"".join(lseq[firstNt:lastNt+1])
                  return (seqCuted,self.findPairGen(secStructCutted))

      def isTerminal(self,i):
          return (self.stableRTab[i+1][0][1] >  self.stableRTab[i][0][1])


      def findPartner(self,n,tab):
          acc = 1
          if(len(tab)>0):
                  stList = list(tab[0].struct)
                  for i in range(n+1,len(stList)):
                          if stList[i] == '(':
                                      acc += 1
                          if stList[i] == ')':
                                      acc -= 1
                                      #print ("index : "+str(i)+" -> "+str(acc))
                          if acc == 0 :
                                      return i
          print ("no match")
          return -1 
        
      def findPartnerFromStruct(self,n,secStruct):
          acc = 1
          stList = list(secStruct)
          for i in range(n+1,len(stList)):
                  if stList[i] == '(':
                          acc += 1
                  if stList[i] == ')':
                          acc -= 1
                          #print ("index : "+str(i)+" -> "+str(acc))
                  if acc == 0 :
                          return i
          print ("no match2")
          print (secStruct + " : "+str(n))
          return -1


      def populateShapeDict(self):
          for ste in self.sESTabsubOpt:
                  if ste.absShape in self.dictShape:
                          self.dictShape[ste.absShape] += 1
          else:
                  self.dictShape[ste.absShape] = 1
          num = 0
          for ste in self.sESTabsubOpt:
                  if not (ste.absShape in self.dictShapeNum):
                          self.dictShapeNum[ste.absShape] = num
                          num += 1
                          #print("test2 : " + str(self.dictShapeNum[ste.absShape]))
                          #else:
          num = 0
          for shape in self.dictShapeNum:
                  if not(shape in self.dictNumShape):
                          self.dictNumShape[self.dictShapeNum[shape]] = ste.absShape
                          #print("test ="+str(self.dictShapeNum[shape]))
                  #print(str(self.dictNumShape[shape]))
                  #else:
          
          
          
      def printShapeDict(self):
          print (json.dumps(self.dictShape, indent=4, sort_keys=True))


      def findPairSo(self):
          if(len(self.sESTabsubOpt)>0):
            self.pairTabSo = []
            self.fillPairTabFirst(self.pairTabSo,self.sESTabsubOpt)

      def findPairMcff(self):
          self.pairTabMcff = []
          self.fillPairTabFirst(self.pairTabMcff,self.sESTabMcff)


      def fillPairTabFirst(self,tab,sESTab):
          if(len(self.sESTabMcff)>0):
                #print self.SEnergyTabMcff[0].struct
                stList = list(sESTab[0].struct)
                #Pour tous les nucleotides (, trouve sa paire
                for n in range(0,len(stList)):
                        if stList[n] == '(':
                                i = self.findPartnerFromStruct(n,sESTab[0].struct)
                                tab.append((n,i)) 
                                
      
      
      #calcul la distance entre la mfe et le tableau de score en le discretisant a l'aide d'un seuiol
      
      def calculDistanceTab(self,mfeTab,probTabtab,seuilTab):
          scroreFin = []
          if(len(mfeTab) == len(probTabtab) and len(mfeTab) == len(seuilTab)):
                  for i in range(0,len(mfeTab)):
                          scroreFin.append(self.calculDistance(mfeTab[i],probTabtab[i],seuilTab[i]))
          return scroreFin

      def calculDistance(self,mfe,probTab,seuil):
          numberOfsiteProbed = 0
          score = 0
          mfeList = list(mfe)
          #make sur that mfe and probTab have the same length
          print ("l_mfe = "+str(len(mfeList))+" ; l_probTab = "+str(len(probTab)))
          for i in range(0,len(probTab[1])):
                  valueSc = 2
                  numberOfsiteProbed +=1
                  #4 case
                  if ((mfe[i] == '(' or mfe[i] == ')') and float(probTab[1][i]) >= seuil):
                          score += 1
                          valueSc = 1
                  if ((mfe[i] == '.') and float(probTab[i]) >= seuil):
                          score -= 1
                          valueSc = -1
                  if ((mfe[i] == '(' or mfe[i] == ')') and float(probTab[1][i]) < seuil):
                          score -= 1
                          valueSc = -1
                  if ((mfe[i] == '.') and float(probTab[1][i]) < seuil):
                          score += 1
                          valueSc = 1 
                  self.distanceScoreTab.append(valueSc)
          finalSc = 0
          if (valueSc == 2):
                  print ("non handled case!!!!!")
          if (float(numberOfsiteProbed) > 0 and float(numberOfsiteProbed)/len(probTab[1])>0.01):
                  finalSc = score / float(numberOfsiteProbed)
                  if (finalSc < 0) :
                          with open(self.negFile, 'a') as f :
                                  f.write(self.name+"\n")
                                  f.write( "prS : " + ";".join(["%.2f" % float(number) for number in probTab[1]]))
                                  f.write( "mfe : " + "".join(mfe)+"\n")
                                  f.write( "seq : " + self.seq+"\n")
                                  f.write( "sco : " + str(finalSc)+"\n")
                                  print ("neg")
          else :
                  with open(self.posFile, 'a') as f :
                          f.write(self.name+"\n")
                          f.write( "prS : " + ";".join(["%.2f" % float(number) for number in probTab[1]]))
                          f.write( "mfe : " + "".join(mfe)+"\n")
                          f.write( "seq : " + self.seq+"\n")
                          f.write( "sco : " + str(finalSc)+"\n")
                          print ("pos")
          if (finalSc == 0):
                  with open(self.zeroFile, 'a') as f :
                          f.write(self.name+"\n")
                          f.write( "prS : " + "|\n".join([",".join(["%.2f" % float(x) for x in probTab[1]])]))
                          f.write( "mfe : " + "".join(mfe)+"\n")
                          f.write( "seq : " + self.seq+"\n")
                          f.write( "sco : " + str(finalSc)+"\n")
                          print ("zero")
              
          else:
                  finalSc = 0
          #print ("scoretab : " + ";".join([str(number) for number in self.distanceScoreTab]))
          #print ("ScoreGlobal : "+str(score ) + "\nmoyenne = "+str(finalSc))
          self.globalScore = finalSc
          return score
      
      def countStableRegion(self,pairNeeded):
          acc = 0
          stablePart = 0
          tempTab = []
          sRFound = False
          for i in range(0,len(self.pairTabSo)-1):
                  if(self.pairTabSo[i][0]+1 == self.pairTabSo[i+1][0] and self.pairTabSo[i][1]-1 == self.pairTabSo[i+1][1]):
                          tempTab.append((self.pairTabSo[i][0],self.pairTabSo[i][1]))
                          acc+=1
                          if (acc == pairNeeded-1):
                                  sRFound = True
                  else:
                          acc = 0
                          if(sRFound):
                                  self.stableRTab.append(tempTab)
                                  stablePart +=1
                          sRFound = False
                          tempTab = []
          if (sRFound):
                  self.stableRTab.append(tempTab)
                  stablePart +=1
          return stablePart

      def createStableRegion(self):
          for regionI in range(0,len(self.stableRTab)):
                  #print(regionI)
                  cuted = self.cutRNA(regionI,self.sESTabsubOpt[0].struct)
                  #print("cuted[0]"+"".join(cuted[0]))
                  #print("cuted[1]"+cuted[1])
                  #RNA2D(self.name+"_"+str(regionI),"".join(cuted[0]),cuted[1],self.fastaPath,"full",False)

      def cutRNA(self,i,secStruct):
          secStructTrans1 = secStruct.replace("[",".")
          secStructTrans = secStructTrans1.replace("]",".")
          firstNt = self.stableRTab[i][0][0]+2
          lastNt = self.stableRTab[i][0][1]
          lseq = list(self.seq)
          #print(secStruct)
          lsecStruct = list(secStructTrans)
          if(len(self.stableRTab) < i+2 or self.isTerminal(i)):
                  lsecStruct[firstNt:lastNt+1]
                  secStructCutted = secStructCutted = "".join(lsecStruct)
                  #print "secStructCutted : "+secStructCutted
                  #print "seqCutted : "+"".join(lseq[firstNt:lastNt+1])
                  #return ("".join(lseq[firstNt:lastNt+1]),self.findPairGen(secStructCutted))
                  return ("".join(lseq[firstNt:lastNt+1]),secStructCutted)
          else:
                  firs2tNt = self.stableRTab[i+1][-1][0]
                  last2Nt = self.stableRTab[i+1][-1][1]
                  secStructCutted = "".join(lsecStruct[firstNt:firs2tNt+1])+"...."+"".join(lsecStruct[last2Nt:lastNt+1])
                  seqCuted = "".join(lseq[firstNt:firs2tNt+1]) +"GAAA"+"".join(lseq[last2Nt:lastNt+1])
                  #print "secStructCutted : "+secStructCutted
                  #print "seqCutted : "+"".join(lseq[firstNt:lastNt+1])
                  #return (seqCuted,self.findPairGen(secStructCutted))
                  return (seqCuted,secStructCutted)

      def isTerminal(self,i):
          return (self.stableRTab[i+1][0][1] >  self.stableRTab[i][0][1])


      def findPartner(self,n,tab):
          acc = 1
          if(len(tab)>0):
                  stList = list(tab[0].struct)
                  for i in range(n+1,len(stList)):
                          if stList[i] == '(':
                                  acc += 1
                          if stList[i] == ')':
                                  acc -= 1
                                  #print ("index : "+str(i)+" -> "+str(acc))
                          if acc == 0 :
                                  return i
          print ("no match")
          return -1 
        
      def findPartnerFromStruct(self,n,secStruct):
          acc = 1
          stList = list(secStruct)
          for i in range(n+1,len(stList)):
                  if stList[i] == '(':
                          acc += 1
                  if stList[i] == ')':
                          acc -= 1
                      #print ("index : "+str(i)+" -> "+str(acc))
                  if acc == 0 :
                          return i
          print ("no match2")
          print (secStruct + " : "+str(n))
          return -1

      def findPartnerGen(self,n,secStruct):
          stList = list(secStruct)
          oppose2 = createDictOpenToClose()
          oppose = createDictCloseToOpen()
          acc = {}
          for key in oppose2:
                  acc[key] = 0
          #En realit/ on a besoin que d'un compteur
          
          stList = list(secStruct)
          acc[stList[n]] = 1
          for i in range(n+1,len(stList)):
                  if stList[i] != '.':
                          if stList[i] == stList[n]:
                                  acc[stList[i]] += 1
                          if stList[i] in oppose :
                                  if stList[n] == oppose[stList[i]]:
                                              acc[oppose[stList[i]]] -= 1
                                              #print ("index : "+str(i)+" -> "+str(acc))
                          if acc[stList[n]] == 0 :
                                  return i
          print ("no match2 : " +stList[n])
          print (secStruct + " : "+str(n))
          return -1
        

        
      def output_files (self,t):
          self.namew = self.name.split('/')[0]
          print ("NAME *****************: " +self.jsFileout+t+".js")
          print ("stat/"+self.specificPath)
          with open(self.jsFileout+t+".js", 'w') as f :
                  f.write(self.jsOut)
                  f.close
          with open(self.htmlFileout+t+".html", 'w') as f :
                  f.write(self.htmlOut)
                  f.close

      def findMotif(self,motif):    
          mo = [m.start() for m in re.finditer(motif,self.seq)]
          iMotif = []
          for m in mo:
                  for k in range(0,len(motif)):
                          iMotif.append(m+k)
          return iMotif



      def createAuxPairTab(self,auxStruct):
          auxPairTabbyPos = []
          auxL = list(auxStruct[0].struct)
          for i in range(0,len(auxStruct[0].struct)):
                  auxPartner = 0
                  if (auxL[i] == "("):
                          auxPartner = self.findPartnerFromStruct(i,auxStruct[0].struct)
                          if (auxPartner != -1):
                                  auxPairTabbyPos.append([auxPartner])
                          else:
                                  print("no partner")
                  else:
                          auxPairTabbyPos.append([])
                          auxL = []
          return auxPairTabbyPos


      
      def jsDictGeneratorByPredictor(self,mcffStr,soStr):
          d_temp = {}
          d_temp["mcff"] = mcffStr
          d_temp["rnaSubOpt"] = soStr
          return d_temp
          
          
  

      

      def createCompoStruct(self,structTab):
          #La structure composite est faite des paires de base les plus observees dans tous les sous-optimaux
          #On a besoin d'un cut-off pour savoir dans combien de structure on voit cette paire de base
          # Dans la 1ere version, le cut-off est a 50% des structures
          
          tabOfPairTab = []
          
          
          
          for struct in structTab:
                  #print("struct : "+ str(struct))
                  tabOfPairTab.append(self.findPairGen(struct.struct))
          
          #print("qqqq len of tabOfPairTab : "+str(len(tabOfPairTab)))

          totalStruct = len(tabOfPairTab)
          dPairCount = {}
          for pairTab in tabOfPairTab:
                  for pair in pairTab:
                          pairStr = "{\"source\":"+str(pair[0])+",\"target\":"+str(pair[1])+",\"value\":1}" #deja formater pour le tableau de lien javascript
                          if(pairStr in dPairCount):
                                  dPairCount[pairStr] = dPairCount[pairStr] + 1
                          else:
                                  dPairCount[pairStr] = 1
                      
          
          
          resultPairTab = []
          for key in dPairCount:
                  #print("key : "+key+" , value: "+str(dPairCount[key]/totalStruct))
                  if (dPairCount[key]/totalStruct > 0.1):
                          resultPairTab.append(key)
                      
                      
          return  resultPairTab
      
     


      def nodeT(self,nom,energy,predictor,p_int):
          d = {}
          d["name"] = nom
          d["predictor"] = predictor
          d["energy"] = energy
          d["p_int"] = p_int
          return d


      def linkT(self,i1,i2,value,pred):
                d = {}
                d["source"] = i1
                d["target"] = i2
                d["value"] = value
                d["pred"] = pred
                #print("linkT")
                return d

      def charTogroup(self,c):
          c = c.upper()
          #print str(c)
          group = ''
          if c == 'A':
                  group = 0
          if c == 'T':
                  group = 1
          if c == 'U':
                  group = 1
          if c == 'C':
                  group = 2
          if c == 'G':
                  group = 3
          return group


      def displayCount(self):
          print ("Nombre de structure sous-optimales de MCff: %d et de RNAsubOpt : %d"  % (len(self.SEnergyTabMcff),len(self.sESTabsubOpt)))

      def restartCount():
          RNA2Daio.count=0

      def displayRNA2D(self):
          print ("Sequence : ", self.seq)
          if len(self.SEnergyTab) != 0  :
                  for s in self.SEnergyTab:
                          #print "type =",type(s)
                          s.displayStructEnergy()
          else :
                  print ("longueur Energy tab nul")
      

  #-------------------------------------------fastaMulti


      def fastaMultiTos(self,filePath):
          inFile  = filePath
          outFile = inFile +".out" 
          print (">> Opening FASTA file...")
          # Reads sequence file list and stores it as a string object. Safely closes file:
          try:  
                  with open(inFile,"r") as newFile:
                          sequences = newFile.read()
                          sequences = re.compile("^>", flags=re.MULTILINE).split(sequences) # Only splits string at the start of a line.
                          del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
                                          # Del removes this empty element.
                          newFile.close()
          except IOError:
                  print ("Failed to open " + inFile)
                  exit(1)

          print (">> Converting FASTA file from multiline to single line and writing to file.")
          # Conversts multiline fasta to single line. Writes new fasta to file.
          for fasta in sequences:
                  try:
                          header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
                  except ValueError:
                          print ("fasta : "+fasta)
                  self.numberOf_windows = int(math.floor(len(sequence)/self.step))
                  for r in range (0,self.numberOf_windows-self.numStep):
                          with open(outFile+str(r)+".fa","w") as newFasta:
                                  print ("lenSeq : "+str(len(sequence))+" r : "+str(r))
                                  header_new = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                                  s = sequence.replace("\n","")
                                  sequence = s.replace('\r','') + "\n" # Replace newlines in sequence, remember to add one to the end.
                                  smallest = len(sequence)
                                  if ((r+self.numStep)*self.step < len(sequence)):
                                          smallest = ((r+self.numStep)*self.step)
                                  self.smallestTab.append(smallest)
                                  print ("r to smallest : ("+str(smallest)+","+str((r+self.numStep)*self.step)+")")
                                  newFasta.write(header_new + sequence[(r*self.step):smallest])
                                  newFasta.close()
                                  print (">> Done! "+outFile+str(r)+".fa")
        
      #processing

      def createFoldConDict(self):
          self.score = []
          if (len(self.SEnergyTabMcff) > 0 and len(self.sESTabsubOpt)>0):
                  mfeSO = self.sESTabsubOpt[0].energy
                  #print "--MFE SO: "+ str(mfeSO)
                  mfeMcff = self.SEnergyTabMcff[0].energy
                  #print "--MFE MCff: "+ str(self.SEnergyTabMcff[0].energy)
                  #donne tous les motifs de 5 possibles dans un dictionnaire
                  fiveNtDi = numberOf_windows(self.seq)
                  #creation des mask  
                  mDict = getMaskMcff(self.seq,fiveNtDi)
                  soDict = getMaskRso(self.seq, fiveNtDi)
          
          for w in fiveNtDi:
                  stri = 'mcff'+str(" -s " + "'"+self.seq+ "'")+str("")+str(" -um "+mDict[w])
                  #print stri + " | "+w
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate() #now you should see your output
                  mcffOutTab = str(output).split("\n")
                  #print "length mcff : "+str(len(mcffOutTab))
                  #print "mcff 1 :"+str(mcffOutTab[0])
                  mfeSplited = mcffOutTab[0].split(' ')
                  condFMcff = float(mfeSplited[1])
                  #print "Mcff MFE ="+ str(condFMcff)
                  difftoperfectMcff = condFMcff - mfeMcff
                  self.condfoldDictMcff[w] = difftoperfectMcff
                  print ("**diffToPerfectMcff : "+str(difftoperfectMcff))
                  
                  maxMcff = 5
                  counter = 0
                  #print "nombre de sous-optimaux: "+str(len(mcffOutTab))
                  for li in mcffOutTab:
                          counter+=1
                          #print li
                          if (counter > maxMcff):
                                  break
                  #print "$$$ " +soDict[w]
                  #creation des fichiers temporaires (mask et seq) 
                  Process=Popen(["echo '"+self.seq+"' > temp"],shell=True,stdout=PIPE,stderr=PIPE)
                  Process=Popen(["echo '"+soDict[w]+"' >> temp"],shell=True,stdout=PIPE,stderr=PIPE)
                  
                  #stri = "RNAsubopt -e 10 −−constraint -s < temp"
                  stri = "RNAsubopt -C -e 1 -s < temp"
                  Process=Popen([stri],shell=True,stdout=PIPE,stderr=PIPE)
                  output,err =  Process.communicate()
                  rsoOutTab = str(output).split("\n")
                  max_so = 5
                  counter=0
                  difftoperfectSO = 0
                  if(len(rsoOutTab)>1):
                          #print "MFE SO :"+ str(rsoOutTab[1])
                          mfeSOsplited = rsoOutTab[1].split()
                          condFSo = float(mfeSOsplited[1])
                          #print "energy : "+mfeSOsplited[1]
                          difftoperfectSO = condFSo - mfeSO
                          print ("**diffToPerfectSO : "+str(difftoperfectSO))
                  for row in rsoOutTab:
                          counter+=1
                          #print "RSo--------> "+row
                          if (counter > max_so):
                                  #print "nombre de sous-optimaux: "+str(len(rsoOutTab))
                                  break
                  #RNAmaskedT.append(RNA2D(self.seq,mTab[counter],"".join(w),mcffOutTab,rsoOutTab))
                  #print "first row : "+ [0]
                  print ("end-----><")
                  self.condfoldDictSO[w] = difftoperfectSO
                  self.score.append(Score(mfeMcff,condFMcff,mfeSO,condFSo,mDict[w],soDict[w],w))
              
          return 0
      
      def correlationFolding(self):
          st=""
          for sc in self.score:
                  st += ","+str(sc.diffMcff)
          st = st[1:]
          f = open("corr_Mcff_"+self.name,'w')
          f.write(st)
          f.close()
          
          st=""
          for sc in self.score:
                  st += ","+str(sc.diffSO)
          st = st[1:]
          f = open("corr_SO_"+self.name,'w')
          f.write(st)
          f.close()
          return 0
      
      
      def commonStruct(self):
          count = 0
          #should pass througth the smallest tab 
          for s in self.sESTabsubOpt:
                  if s.struct in self.mcffD:
                          count = count + 1
          return count 


      def adMat2(self,sESTab):
          # 'd' contiendra un dictionnaire de dictionnaire 
          d = {}
          #pour toute les sous-optimaux, on fait la liste des structure voisine (exponentiel ... sur la longueur de la sequence)
          for so in sESTab:
                      soPairTab1 = self.findPairGen( so.struct )
                      lPairTab = len(soPairTab1)
                      if(lPairTab >0):
                              d[so.struct] = {}
                              for so2 in sESTab:
                                  soPairTab2 = self.findPairGen( so2.struct )
                                  lPairTab = max(len(soPairTab2),lPairTab)
                                  ressemblance = self.commonPair(soPairTab1,soPairTab2)/float(lPairTab)
                                  #print ("resemblance = "+str(ressemblance))
                                  if(ressemblance > 0.95):
                                          #if(so.struct in self.indexTabMcff and so2.struct in self.indexTabMcff):
                                              #print (str(self.indexTabMcff[so.struct])+","+str(so2.energy - so.energy)+","+str(self.indexTabMcff[so2.struct]))
                                          d[so.struct][so2.struct]     =    ressemblance
                        #print "Total : "+so.struct +" : "+ str(len(dtemp)) + " sur : "+ str(len(sESTab))
          return d

      def commonPair(self,so1Tab,so2Tab):
          count = 0
          for p1 in so1Tab:
                      for p2 in so2Tab:
                              if(p1[0] == p2[0] and p1[1] == p2[1]):
                                  #print ("common pair :"+str(p1[0])+","+str(p2[0])+","+str(p1[1])+","+str(p2[1]))
                                  count += 1
          return count

      def createAdList(self,indexTab,adjMat):
          tab = []
          for s in adjMat:
                      #print "-----"+name
                      for s2 in adjMat[s]:
                              if(s != s2):
                                #print (str(str(indexTab[s])+";"+str(adjMat[s][s2])+";"+str(indexTab[s2])))
                                if (-1 *adjMat[s][s2]< 0):
                                    #f.write(s+","+str(adjMat[s][s2])+","+s2+"\n")
                                    tab.append((indexTab[s],adjMat[s][s2],indexTab[s2]))
          #f.close()
          return tab

        



      def createAdjListSO(self,path):
          self.adjListSO = []
          if not os.path.exists(path+"/adList/"):
                      os.makedirs(path+"/adList/")
         
          for s in self.adjMatSO:
                      #print "-----"
                      for s2 in self.adjMatSO[s]:
                              if(s != s2):
                                first = str(self.indexTabSO[s])
                                second = str(self.adjMatSO[s][s2])
                                third = str(self.indexTabSO[s2])
                                #print ("createAdjListSO : "+first+";"+second+";"+third)
                                if (-1 * self.adjMatSO[s][s2]<0):
                                    if(s in self.indexTabSO and s2 in self.indexTabSO):
                                            #f.write(s+","+str(self.adjMatSO[s][s2])+","+s2+"\n")
                                            self.adjListSO.append((self.indexTabSO[s],self.adjMatSO[s][s2],self.indexTabSO[s2]))
       
    
      def printDict(self):
          sorted_x = sorted(self.condfoldDictMcff.items(), key=operator.itemgetter(1))
          f = open('mcff.csv', 'w') 
          for tu in sorted_x:
                      f.write(tu[0]+","+str(tu[1])+"\n")
          f.close()
          sorted_y = sorted(self.condfoldDictSO.items(), key=operator.itemgetter(1))
          f = open('SO.csv', 'w') 
          for tu in sorted_y:
                      f.write(tu[0]+","+str(tu[1])+"\n")
          f.close()
      
      #def composanteConnexes(self):
        #tabFinale = []
        ##parcours de toutes les structures
        #for s in self.adjMat:
              #break
          ##creation d'un noeud 



      #mAdToD3Graph prend en entree deux listes d'adjencences et deux listes de conformations et donne un graph de transition en format noeud + link
      def mAdToD3Graph(self,sEtabMCff,sEtabSO,adjListMcff,adjListSO):
          maxMcff = max([x.energy for x in sEtabMCff])
          minMcff = min([x.energy for x in sEtabMCff])
          maxSO = max([x.energy for x in sEtabSO])
          minSO = min([x.energy for x in sEtabSO])
          nodeTempsT = ""
          tabNode = []
          for i in range(0,len(sEtabMCff)):
                      tabNode.append(self.nodeT(i,sEtabMCff[i].energy,"mcff",1))

          for i in range(0,len(sEtabSO)):
                      tabNode.append(self.nodeT(i+len(sEtabMCff),sEtabSO[i].energy,"so",2))

          linkTempsT = ""
          tabLink = []
          print("adjList len : "+str(len(sEtabMCff)))
          for triplet in adjListMcff:
                      tabLink.append( self.linkT(triplet[0],triplet[2],triplet[1],"mcff"))#(source,target,value)
                              
          for triplet in adjListSO:
                      tabLink.append( self.linkT(triplet[0]+len(sEtabMCff),triplet[2]+len(sEtabMCff),triplet[1],"so"))#(source,target,value)

          for k in self.sESTabMcff:
                      s1 = k.struct
                      for k2 in  self.sESTabsubOpt:
                                s2 = k2.struct
                                #print("s1:"+s1+"\ns2:"+s2)
                                soPairTab2 = self.findPairGen( s2 )
                                soPairTab1 = self.findPairGen( s1 )
                                lPairTab = max(len(soPairTab2),len(soPairTab1))
                                ressemblance = self.commonPair(soPairTab1,soPairTab2)/float(lPairTab)
                                if(ressemblance > 0.8):
                                          print("2 structures semblable!")
                                          #print("so : "+str(self.indexTabSO[s2]))
                                          #print("mcff : "+str(self.indexTabMcff[s1]))
                                          tabLink.append( self.linkT(self.indexTabSO[s2]+len(sEtabMCff),self.indexTabMcff[s1],1,"so_mcff"))#(source,target,value)

          
          #print ("tabLink : "+str(tabLink[:100]))
          return (tabNode,tabLink,maxMcff,minMcff,maxSO,minSO)






      def displayCount(self):
          print ("Nombre de structure sous-optimales de MCff: %d et de RNAsubOpt : %d"  % (len(self.SEnergyTabMcff),len(self.sESTabsubOpt)))

      def restartCount():
          RNA2Daio.count=0

      def displayRNA2D(self):
          print ("Sequence : ", self.seq)
          if len(self.SEnergyTab) != 0  :
                  for s in self.SEnergyTab:
                          #print "type =",type(s)
                          s.displayStructEnergy()
          else :
                  print ("longueur Energy tab nul")
      

      #-------------------------------------------fastaMulti


      def fastaMultiTos(self,filePath):
          inFile  = filePath
          outFile = inFile +".out" 

          print (">> Opening FASTA file...")
          # Reads sequence file list and stores it as a string object. Safely closes file:
          try:  
                  with open(inFile,"r") as newFile:
                          sequences = newFile.read()
                          sequences = re.compile("^>", flags=re.MULTILINE).split(sequences) # Only splits string at the start of a line.
                          del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
                                  # Del removes this empty element.
                  newFile.close()
          except IOError:
                  print ("Failed to open " + inFile)
                  exit(1)

          print (">> Converting FASTA file from multiline to single line and writing to file.")
          # Conversts multiline fasta to single line. Writes new fasta to file.
          for fasta in sequences:
                  try:
                              header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
                  except ValueError:
                              print (fasta)
                  self.numberOf_windows = int(math.floor(len(sequence)/self.step))
                  for r in range (0,self.numberOf_windows-self.numStep):
                              with open(outFile+str(r)+".fa","w") as newFasta:
                                      print ("lenSeq : "+str(len(sequence))+" r : "+str(r))
                                      header_new = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                                      s = sequence.replace("\n","")
                                      sequence = s.replace('\r','') + "\n" # Replace newlines in sequence, remember to add one to the end.
                                      smallest = len(sequence)
                                      if ((r+self.numStep)*self.step < len(sequence)):
                                              smallest = ((r+self.numStep)*self.step)
                                      self.smallestTab.append(smallest)
                                      print ("r to smallest : ("+str(smallest)+","+str((r+self.numStep)*self.step)+")")
                                      newFasta.write(header_new + sequence[(r*self.step):smallest])
                                      newFasta.close()
                                      print (">> Done! "+outFile+str(r)+".fa")
      
      def createTabD_nt(self,tabNt,scoreLabel):
          tabResult = []
          id_u = 0
          for nt in range(0,len(tabNt)):
            
              id_u +=1
              tabResult.append({})
              tabResult[nt]["cpjp"] = 0
              tabResult[nt]["ncm"] = {}
              for ncm in tabNt[nt]["localNcmD"]['mcff']:
                      tabResult[nt]["ncm"][ncm+"_mcff"] = tabNt[nt]["localNcmD"]['mcff'][ncm]
              #print('len(tabNt[nt]["localNcmD"]["mcff"]) :'+str(len(tabNt[nt]["localNcmD"]['mcff'])))
              
              for ncm in tabNt[nt]["localNcmD"]['so']:
                      tabResult[nt]["ncm"][ncm+"_so"] = tabNt[nt]["localNcmD"]['so'][ncm]
              
                
            #print('len(tabNt[nt]["localNcmD"]["so"]) :'+str(len(tabNt[nt]["localNcmD"]['so'])))
            #print(str(tabNt[nt]["sorte"]))
              tabResult[nt]["sorte"]  =  tabNt[nt]["sorte"] 
              tabResult[nt]["score"]  =  tabNt[nt]["score"]
              tabResult[nt]["erreur"]  =  tabNt[nt]["erreur"]
              tabResult[nt]["rna_id"] = tabNt[nt]["rna_id"]
              tabResult[nt]["position"] = tabNt[nt]["position"]
              tabResult[nt]["inHelice"] = tabNt[nt]["inHelice"]
              tabResult[nt]["freqPairee_mcff"] = tabNt[nt]["freqPairee_mcff"]
              tabResult[nt]["freqPairee_so"] = tabNt[nt]["freqPairee_so"]
              tabResult[nt]["scoreLabel"] = scoreLabel
              tabResult[nt]["u_id"] = [id_u]
          #print("tabResult : "+str(len(tabResult)))
          return tabResult
            
      #-------------------------------------------------------------------------      
      def basicPath(self):
          d = {}
           
          d["score"]  =  ["nts","i","score"]
          d["erreur"]  =  ["nts","i","erreur"]
          d["rna_id"] = ["rna_id"]
          d["position"] = ["position"]
          
          return d
      #-------------------------------------------------------------------------  
      def sorteHeaderGen(self):
          t = ["u_id","rna_id","position","score","erreur","sorte","scoreLabel"]
          return t
      
      def sortePathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          return d
        
      def dictForSorte(self):
          name = "sorte"
          header = self.sorteHeaderGen()
          path1 = self.sortePathGen()
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d
          
          
          
          
      def seqMotif3HeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","seqMotif3"]
          return t
      
      def seqMotif3PathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["seqMotif3"] = ["nts","i","seqMotif3"]
          return d
      
      def dictForSeqMotif3(self):
          name = "seqMotif3"
          header = self.seqMotif3HeaderGen()
          path1 = self.seqMotif3PathGen()
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d
        
        
        
      def seqMotif5HeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","seqMotif5"]
          return t
      
      def seqMotif5PathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["seqMotif5"] = ["nts","i","seqMotif5"]
          return d
      
      def dictForSeqMotif5(self):
          name = "seqMotif5"
          header = self.seqMotif5HeaderGen()
          path1 = self.seqMotif5PathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          d  = {"name":name,"header":header,"path":path2}
          return d      
      
      
      
      def basicHeaderGen(self):
          t = ["u_id","score","erreur","scoreLabel","rna_id","position"]
          return t
      
      def basicPathGen(self):
          tab = []
          d = self.basicPath()
          return d

      def dictForBasic(self):
          name = "basic"
          header = self.basicHeaderGen()
          path1 = self.basicPathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d 
      
      def pairPartnerHeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","mfp_mcff","mfp_so"]
          return t
          
      def pairPartnerPathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["mfp_mcff"] = ["nts","i","mfp_mcff"]
          d["mfp_so"] = ["nts","i","mfp_so"]
          return d

      def dictForPairPartnerPathGen(self):
          name = "pairPartner"
          header = self.pairPartnerHeaderGen()
          path1 = self.pairPartnerPathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d 
      
      
      
      
      
      def regionHeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","region"]
          return t    
      
      def regionPathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["region"] = ["nts","i","region"]
          return d

      def dictForRegion(self):
          name = "region"
          header = self.regionHeaderGen()
          path1 = self.regionPathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d 
        
        
        
       
      
      def toToWL(self,d):
        return " ".join(d.keys())
      
      
      def mfeNcmHeaderGen_so(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_mfe_so"] 
          return t    
      
      def mfeNcmPathGen_so(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_mfe_so"] = ["nts","i","localNcmD_mfe","so"]
          return d

      def dictForMfeNcm_so(self):
          name = "mfeNcm_so"
          header = self.mfeNcmHeaderGen_so()
          path1 = self.mfeNcmPathGen_so()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['so_ncm']
          d  = {"name":name,"header":header,"path":path2,"listC":"so_ncm"}
          return d 
          
          

      
      def mfeNcmHeaderGen_mcff(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_mfe_mcff"]
          return t    
      
      def mfeNcmPathGen_mcff(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_mfe_mcff"] = ["nts","i","localNcmD_mfe","mcff"]
          return d

      def dictForMfeNcm_mcff(self):
          name = "mfeNcm_mcff"
          header = self.mfeNcmHeaderGen_mcff()
          path1 = self.mfeNcmPathGen_mcff()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1]  + self.listC['mcff_ncm']
          d  = {"name":name,"header":header,"path":path2,"listC":"mcff_ncm"}
          return d           
          
          

        
        
      def ncm25HeaderGen_so(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_so"] 
          return t    
      
      def ncm25PathGen_so(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_so"] = ["nts","i","localNcmD","so"]
          return d
        
      def dictForNcm25_so(self):
          name = "ncm25_so"
          header = self.ncm25HeaderGen_so()
          path1 = self.ncm25PathGen_so()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['so_ncm']
          d  = {"name":name,"header":header,"path":path2,"listC":"so_ncm"}
          return d 
        
          
          
      def ncm25HeaderGen_mcff(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_mcff"]
          return t    
      
      def ncm25PathGen_mcff(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_mcff"] = ["nts","i","localNcmD","mcff"]
          return d
        
      def dictForNcm25_mcff(self):
          name = "ncm25_mcff"
          header = self.ncm25HeaderGen_mcff() 
          path1 = self.ncm25PathGen_mcff()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['mcff_ncm']
          d  = {"name":name,"header":header,"path":path2,"listC":"mcff_ncm"}
          return d                 
        
        
      def ncm_mfe_detailHeaderGen_so(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_mfe_detail_so"] 
          return t    
      
      def ncm_mfe_detailPathGen_so(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_mfe_detail_so"] = ["nts","i","localNcmD_mfe_detail","so"]
          return d
      

      def dictForNcm_mfe_detail_so(self):
          name = "ncm_mfe_detail_so"
          header = self.ncm_mfe_detailHeaderGen_so() 
          path1 = self.ncm_mfe_detailPathGen_so()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['so_detail']
          d  = {"name":name,"header":header,"path":path2,"listC":"so_detail"}
          return d 

        
        
      def ncm_mfe_detailHeaderGen_mcff(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_mfe_detail_mcff"] 
          return t    
      
      def ncm_mfe_detailPathGen_mcff(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_mfe_detail_mcff"] = ["nts","i","localNcmD_mfe_detail","mcff"]
          return d
      

      def dictForNcm_mfe_detail_mcff(self):
          name = "ncm_mfe_detail_mcff"
          header = self.ncm_mfe_detailHeaderGen_mcff()
          path1 = self.ncm_mfe_detailPathGen_mcff()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['mcff_detail']
          d  = {"name":name,"header":header,"path":path2,"listC":"mcff_detail"}
          return d 

    


      def ncm_detailHeaderGen_so(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_detail_so"] 
          return t    
      
      def ncm_detailPathGen_so(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_detail_so"] = ["nts","i","localNcmD_detail","so"]
          return d
      

      def dictForNcm_detail_so(self):
          name = "ncm_detail_so"
          header = self.ncm_detailHeaderGen_so()
          path1 = self.ncm_detailPathGen_so()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['so_detail']
          d  = {"name":name,"header":header,"path":path2,"listC":"so_detail"}
          return d 
        
        
        
        
      def ncm_detailHeaderGen_mcff(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","ncm_detail_mcff"] 
          return t    
      
      def ncm_detailPathGen_mcff(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["ncm_detail_mcff"] = ["nts","i","localNcmD_detail","mcff"]
          return d
      

      def dictForNcm_detail_mcff(self):
          name = "ncm_detail_mcff"
          header = self.ncm_detailHeaderGen_mcff()
          path1 = self.ncm_detailPathGen_mcff()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          header = header[:-1] + self.listC['mcff_detail']
          d  = {"name":name,"header":header,"path":path2,"listC":"mcff_detail"}
          return d         

        
        
        
      '''
        
      def ncmPosHeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","freqPairee_mcff","freqPairee_so"]
          return t    
              
      def ncmPosPathGen(self):
          tab = []
          d = self.basicPath()
          return d

      def dictForSeqMotif5(self):
          name = "basic"
          header = self.basicHeaderGen()
          path1 = self.basicPathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d 
      ''' 
      def ncmVoisinHeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","voisinDroit_sorte","voisinDroit_etat_mcff","voisinDroit_etat_so","voisinGauche_sorte","voisinGauche_etat_mcff","voisinGauche_etat_so"]
          return t   
      
      def ncmVoisinPathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["voisinDroit_sorte"] = ["nts","i","voisinDroit","sorte"]
          d["voisinDroit_etat_mcff"] = ["nts","i","voisinDroit","etat","mcff"]
          d["voisinDroit_etat_so"] = ["nts","i","voisinDroit","etat","rnaSubOpt"]
          d["voisinGauche_sorte"] = ["nts","i","voisinGauche","sorte"]
          d["voisinGauche_etat_mcff"] = ["nts","i","voisinGauche","etat","mcff"]
          d["voisinGauche_etat_so"] = ["nts","i","voisinGauche","etat","rnaSubOpt"]
          return d

      def dictForNcmVoisin(self):
          name = "ncmVoisin"
          header = self.ncmVoisinHeaderGen()
          path1 = self.ncmVoisinPathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d 
      
      def heliceHeaderGen(self):
          t = ["u_id","score","erreur","sorte","scoreLabel","rna_id","position","inHelice"]
          return t   
      
      def helicePathGen(self):
          tab = []
          d = self.basicPath()
          d["sorte"]  =  ["nts","i","sorte"]
          d["inHelice"] = ["nts","i","inHelice"]
          return d

      def dictForHelice(self):
          name = "helice"
          header = self.heliceHeaderGen()
          path1 = self.helicePathGen()
          
          path2 = []
          #To have the same order as header
          for k in header:
            if(k == "u_id"):
              path2.append("u_id")
            else:
              if(k == "scoreLabel"):
                path2.append("scoreLabel")
              else:
                path2.append(path1[k])
          
          
          d  = {"name":name,"header":header,"path":path2}
          return d 
        
        
        
 
          
     
      


      def createCaracteristicHi(self):
          d = {}
          '''
          d["2_5_so"] = 1
          d["2_5_mcff"] = 1
          d["2_6_so"] = 1
          d["2_6_mcff"] = 1
          d["3_5_mcff"] = 1
          d["6_4_so"] = 1
          d["6_so"] = 1
          d["Premier_Paired_mcff"] = 1
          d["Premier_Paired_so"] = 1
          d["Premier_notPaired_mcff"] = 1
          d["Premier_notPaired_so"] = 1
          '''
          return d

      def createCaracteristicLow(self):
          d = {}
          '''
          d["2_5_so"] = 1
          d["2_5_mcff"] = 1
          d["2_6_so"] = 1
          d["2_6_mcff"] = 1
          d["3_5_mcff"] = 1
          d["6_4_so"] = 1
          d["6_so"] = 1
          d["Premier_Paired_mcff"] = 1
          d["Premier_Paired_so"] = 1
          d["Premier_notPaired_mcff"] = 1
          d["Premier_notPaired_so"] = 1
          '''
          return d

      def ncmCaracteristicsInit(self):
          d = {}
          d["cpjp"] = 1
          d["10_so"] = 1
          d["L_so"] = 1
          d["2_2_mcff"] = 1
          d["2_2_so"] = 1
          d["2_3_mcff"] = 1
          d["2_3_so"] = 1
          d["2_4_mcff"] = 1
          d["2_4_so"] = 1
          d["2_5_mcff"] = 1
          d["2_5_so"] = 1
          d["2_6_mcff"] = 1
          d["2_6_so"] = 1
          d["2_L_mcff"] = 1
          d["2_L_so"] = 1
          d["21_so"] = 1
          d["25_so"] = 1
          d["3_2_mcff"] = 1
          d["3_2_so"] = 1
          d["3_3_mcff"] = 1
          d["3_3_so"] = 1
          d["3_4_mcff"] = 1
          d["3_4_so"] = 1
          d["3_5_mcff"] = 1
          d["3_5_so"] = 1
          d["3_6_so"] = 1
          d["3_L_mcff"] = 1
          d["3_L_so"] = 1
          d["4_2_mcff"] = 1
          d["4_2_so"] = 1
          d["4_3_mcff"] = 1
          d["4_3_so"] = 1
          d["4_4_mcff"] = 1
          d["4_4_so"] = 1
          d["4_5_so"] = 1
          d["4_6_so"] = 1
          d["4_L_mcff"] = 1
          d["4_L_so"] = 1
          d["4_mcff"] = 1
          d["5_2_mcff"] = 1
          d["5_2_so"] = 1
          d["5_3_mcff"] = 1
          d["5_3_so"] = 1
          d["5_4_so"] = 1
          d["5_5_so"] = 1
          d["5_6_so"] = 1
          d["5_L_mcff"] = 1
          d["5_L_so"] = 1
          d["5_mcff"] = 1
          d["5_so"] = 1
          d["6_2_mcff"] = 1
          d["6_2_so"] = 1
          d["6_3_so"] = 1
          d["6_4_so"] = 1
          d["6_5_so"] = 1
          d["6_6_so"] = 1
          d["6_L_mcff"] = 1
          d["6_L_so"] = 1
          d["6_mcff"] = 1
          d["6_so"] = 1
          d["7_mcff"] = 1
          d["7_so"] = 1
          d["8_mcff"] = 1
          d["8_so"] = 1
          d["L_2_mcff"] = 1
          d["L_4_so"] = 1
          d["L_5_so"] = 1
          d["L_6_so"] = 1
          d["L_L_mcff"] = 1
          d["L_L_so"] = 1
          d["Premier_Paired_mcff"] = 1
          d["Premier_Paired_so"] = 1
          d["Premier_notPaired_mcff"] = 1
          d["Premier_notPaired_so"] = 1
          d["noVD_Paired_mcff"] = 1
          d["noVD_Paired_so"] = 1
          d["noVD_notPaired_mcff"] = 1
          d["noVD_notPaired_so"] = 1
          d["noVG_Paired_mcff"] = 1
          d["noVG_Paired_so"] = 1
          d["noVG_notPaired_mcff"] = 1
          d["noVG_notPaired_so"] = 1
          d2 = {}
          for k in d.keys():
              if(k not in self.cpjp):
                  d2[k] = 0
          return d2
          



def arrayToDict(a):
  r  = {}
  for e in a:
    r[e] = 1;
  return r

          
def normalise(d,kind):
    if(kind == "sum"):
      s = sum(d.values())
      if(s == 0):
          #print("sum = 0")
          #for i in d.itervalues():
            #print(str(i))
          factor = 0
      else:
          factor=1.0/s
      for k in d:
          d[k] = d[k]*factor
    return d

def createDictOpenToClose():
    oppose2 = {}
    oppose2['(']=')'
    oppose2['{']='}'
    oppose2['[']=']'
    oppose2['<']='>'
    oppose2['A']='a'
    oppose2['B']='b'
    return oppose2 

def createDictCloseToOpen():
    oppose = {}
    oppose[')']='('
    oppose['}']='{'
    oppose[']']='['
    oppose['>']='<'
    oppose['a']='A'
    oppose['b']='B'
    return oppose


def tab2Dict(tab):
  d = {}
  for e in tab:
      d[e.struct]=e.rank
  return d

def printPair(ptab):
  for p in ptab:
          print ("("+str(p[0])+","+str(p[1])+")")

def ntToCsv(fileName,d):
    #my_dict = {"test": 1, "testing": 2}

    with open(fileName, 'w') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, d.keys())
        w.writeheader()
        w.writerow(d)
    return 1

def extend3(d_temp,src):
  newD = {}
  for k in src:
    #print(k)
      if(k in d_temp):
          #print("d_temp[k] : "+str(d_temp[k]))
          #print("src[k] : "+str(src[k]))
          newD[k] = float(d_temp[k]) + src[k]
      else:
          newD[k] = src[k]
          
  return newD








  
def extend_mcff(d,src):
  newD = {}
  for k in src:
      #print(k)
      if(k == 'cpjp'):
          if(k in d ):
              #print("d[k] : "+str(d[k]))
              #print("src[k] : "+str(src[k]))
              newD[k] = float(d[k]) + src[k]
          else:
              newD[k] = src[k]
          
      else:
          if(k[:-4] in d ):
          #print("d[k] : "+str(d[k]))
          #print("src[k] : "+str(src[k]))
              newD[k] = float(d[k[:-4]]) + src[k]
          else:
              newD[k] = src[k]
  return newD

def extend_so(d,src):
  newD = {}
  for k in src:
    #print(k)
      if(k != 'scoreLabel'):
          if(k == 'cpjp'):
              if(k in d ):
                  #print("d[k] : "+str(d[k]))
                  #print("src[k] : "+str(src[k]))
                  newD[k] = float(d[k]) + src[k]
              else:
                  newD[k] = src[k]
              
          else:
              
              if(k[:-3] in d ):
                  #print("d[k] : "+str(d[k]))
                  #print("src[k] : "+str(src[k]))
                  
                  newD[k] = float(d[k[:-3]]) + src[k]
              else:
                  newD[k] = src[k]
  return newD
      
  return newD
def extractKeys(dTab):
    d = {}
    for d_temp in dTab:
      for k in d_temp:
          if(not k in d):
              d[k] = 0
    return d


    
def createNtCsvHeader(filename,keyCsvHeaderTab):
  f = open(filename,'w')
  fK_str = ",".join(keyCsvHeaderTab)
  f.write(fK_str+keys)
  f.close()
#-----------------------------------------------------------------------------------------------------------------------------------------------StructureEnergyTuple
class StructEnergy:
    'structure/Energy'
    count = 0

    def __init__(self, struct,seq, energy, summary,rank):
          self.seq = seq
          self.struct = struct
          self.pairTab = []
          self.findPair()
          #print ("nombre de paire : "+ str(len(self.pairTab)))
          self.tripletTabForJS = []
          self.heliceTab = []
          self.fillHeliceTab()
          #print ("HeliceTab : "+str(self.heliceTab))
          self.heliceTabObj = []
          #print "nombre d'helices : "+ str(len(self.heliceTab))
          self.jonctionTab = []
          self.root = []
          #self.printjunction()
          self.absShape = ""
          '''
          self.createTree()
          #                                                                                                                                                                                                                                                                                                                                                                                      print ("nombre d'helices root : "+ str(len(self.root)))
          for i in range(0,len(self.heliceTabObj)):
                  self.heliceTabObj[i].createJ(None)
          #"\"squelette\":
          tabTemp = []
          self.stringForJSroot = "[" 
          for h in self.root:
                  tabTemp.append( h.toStringForJS())
                  self.absShape += h.shape()
          self.stringForJSroot += ",".join(tabTemp)
          self.stringForJSroot += "]"
          '''
          self.energy = energy
          self.summary =  summary
          self.rank =  rank
          StructEnergy.count += 1
          
          
    def __str__(self):
          return "\nSeq  : "+self.seq+"\nSSec : "+self.struct + "\nEnergy : "+ self.energy
          
    def toString(self):
          return self.struct + " "+ str(self.energy)
    
            
            
    def displayCount(self):
          print ("Nombre de structure 2D %d" % StructEnergy.count)
  
    def restartCount(self):
          StructEnergy.count=0

    def displayStructEnergy(self):
          s = "".join(["Structure: ", self.struct,  ", Energy: ", str(self.energy), " Summary: ", self.summary, " Rank: ", str(self.rank)])
          return s
  
  #processing

    def nearStruct(self):
          adstructTab = self.nearstructHelper(self.struct)
          ad2strucDict = {}
          ad2strucDict.clear()
          for st1 in adstructTab:
                  ad2strucDict.update(self.nearstructHelper(st1))
          ad2strucDictWellP = {}
          ad2strucDictWellP.clear()
          for st2 in ad2strucDict:
                  if (self.struct != st2 and balanced(st2)):
                          ad2strucDictWellP[st2]=0
                          #print st2
          return ad2strucDictWellP


    def nearstructHelper(self,st):
          adstructTab = {}
          stReel = st
          for i in range(0,len(st)):
                  if (stReel[i] == '.'):
                          st = '%s' % stReel
                          stList = list(st)
                          stList[i] = '('
                          newS1 = "".join(stList)
                          #print(newS1)
                          adstructTab[newS1]=1
                          st = '%s' % stReel
                          stList = list(st)
                          stList[i] = ')'
                          newS2 = "".join(stList)
                          adstructTab[newS2]=1
                  if (stReel[i] == ')'):
                          st = '%s' % stReel
                          stList = list(st)
                          stList[i] = '.'
                          newS1 = "".join(stList)
                          adstructTab[newS1]=1
                          st = '%s' % stReel
                          stList = list(st)
                          stList[i] = '('
                          newS2 = "".join(stList)
                          adstructTab[newS2]=1
                  if (stReel[i] == '('):
                          st = '%s' % stReel
                          stList = list(st)
                          stList[i] = '.'
                          newS1 = "".join(stList)
                          adstructTab[newS1]=1
                          st = '%s' % stReel
                          stList = list(st)
                          stList[i] = ')'
                          newS2 = "".join(stList)
                          adstructTab[newS2]=1
          return adstructTab

    def findPartnerFromStruct(self,n):
          acc = 1
          stList = list(self.struct)
          for i in range(n+1,len(stList)):
                  if stList[i] == '(':
                          acc += 1
                  if stList[i] == ')':
                          acc -= 1
                          #print ("index : "+str(i)+" -> "+str(acc))
                  if acc == 0 :
                          return i
          print ("no match")
          print (self.struct + " : "+str(n))
          return -1

    def findPair(self):
          #print "startFindPair"
          #print self.sESTabsubOpt[0].struct
          stList = list(self.struct)
          #Pour tous les nucleotides (, trouve sa paire
          for n in range(0,len(stList)):
                  if stList[n] == '(':
                          i = self.findPartnerFromStruct(n)
                          if(i == -1):
                                  print ("pairNOTfoundSO")
                          self.pairTab.append((n,i))
          #print "pairfound"

    def fillHeliceTab(self):
          sRFound = False
          tempTab = []
          acc = 0
          pairNeeded = 2
          for i in range(0,len(self.pairTab)):
                  if(i < len(self.pairTab)-1 and self.pairTab[i][0]+1 == self.pairTab[i+1][0] and self.pairTab[i][1]-1 == self.pairTab[i+1][1]):
                          tempTab.append((self.pairTab[i][0],self.pairTab[i][1]))
                          acc+=1
                          if (acc == pairNeeded-1):
                                  sRFound = True
                  else:
                          acc = 0
                          if(sRFound):
                                  tempTab.append((self.pairTab[i][0],self.pairTab[i][1]))
                                  self.heliceTab.append(tempTab)
                          sRFound = False
                          tempTab = []
          if (sRFound):
                  self.heliceTab.append(tempTab)




    def createTree(self):
          self.root = []
          self.heliTab = []
          #print("createTree")
          heliTab = []
          #print "self.heliceTab[0] : "+ str(self.heliceTab[0])
          h_courante = Helice(self.heliceTab[0][0][0],self.heliceTab[0][-1][0],self.heliceTab[0][-1][1],self.heliceTab[0][0][1],"",self.heliceTab[0],0)
          self.heliceTabObj.append(h_courante)
          h_T = [h_courante]
          self.root.append(h_courante)
          h_restanteTab = self.heliceTab[1:]
          indice = self.createTreeHelper(h_T,h_restanteTab,1)
          #print("finCreateTree")
          return 1

    def h1_In_h2(self,h2,h1):
          return h2.milieu1I < h1.debutI and h1.finI < h2.milieu2I
        
    #Verifie si la premiere helice restante est dans h_courante sinon pop tabtree et reapelle createTreeHelper
    def createTreeHelper(self,h_T,h_restanteTab,i):
          if len(h_restanteTab) == 0:
                  return
          new_H = Helice(h_restanteTab[0][0][0],h_restanteTab[0][-1][0],h_restanteTab[0][-1][1],h_restanteTab[0][0][1],"",h_restanteTab[0],i)
          self.heliceTabObj.append(new_H)
          #print "Longueur de h_T : "+str(len(h_T))
          l = len(h_T)
          if (l>0):
                  if (self.h1_In_h2(h_T[-1],new_H)):
                          #print "oui"
                          h_T[-1].helices.append(new_H)
                          #print(heli.toString(l))
                          h_T.append(new_H)
                          self.createTreeHelper(h_T,h_restanteTab[1:],i+1)
                  else:
                          #print "non"
                          h_T.pop()
                          self.createTreeHelper(h_T,h_restanteTab,i)
          else:
                  #print("Root"+" : ")
                  #print(heli.toString(l))
                  self.root.append(new_H)
                  h_T.append(new_H)
                  self.createTreeHelper(h_T,h_restanteTab[1:],i+1)

iparens = iter('()')
parens = dict(zip(iparens, iparens))
closing = parens.values()

def balanced(astr):
    stack = []
    for c in astr:
        d = parens.get(c, None)
        if d:
            stack.append(d)
        elif c in closing:
            if not stack or c != stack.pop():
                return False
    return not stack

def test():
  rna2d = RNA2Daio("gggggccccc")
  rna2d.displayRNA2D()

def window(fseq):
  w_size=5
  diAns = {}
  l = list(fseq)
  for i in range(0,len(l)-w_size+1):
          tempTab = []
          for i2 in range(i,i+w_size):
                  tempTab.append(l[i2])
          w = "".join(tempTab)
          if(w in diAns):
                  diAns[w]+=1
          else:
                  diAns[w]=1
  return diAns

def getMaskMcff(seq,di):
  di_ans={}
  w_size=5
  for w in di:
          st = list('x'*len(seq))
          indexTab = [m.start() for m in re.finditer(w,seq)]
          #print "indexTab"+str(indexTab)
          for i in indexTab:
                  for i2 in range(i,(i+w_size)):
                          st[i2] = '.'
          di_ans[w]="".join(st)
  return di_ans

def getMaskRso(seq,di,w_size=5):
  di_ans = {}
  for w in di:
          st = list('.'*len(seq))
          indexTab = [m.start() for m in re.finditer(w,seq)]
          for i in indexTab:
                  for i2 in range(i,i+w_size):
                            st[i2] = 'x'
          di_ans[w]="".join(st)
  return di_ans
  
#----------------------------------------------------------------------------------------------structureEnergy from multibranche3D

class Helice:
    count = 0
    def __init__(self, debutI,milieu1I,milieu2I,finI,seqTotal,pairTab,indice):
          self.indice = indice
          self.pairTab = pairTab
          self.debutI = debutI
          self.milieu1I = milieu1I
          self.milieu2I = milieu2I
          self.finI = finI
          self.seq = seqTotal[debutI:milieu1I+1] + seqTotal[milieu2I:finI+1]
          self.helices = []
          self.index = debutI * len(seqTotal) + finI
          self.junction1 = None
          self.junction2 = None
  
    def toString(self,n):
          stri = str(self.debutI) +"-"+str(self.milieu1I)+"-"+str(self.milieu2I) +"-"+str(self.finI)
          n += 1
          #print(str(n))
          for h in self.helices:
                  stri += "\n"+self.addSpace(n) + h.toString(n)
          return stri
    
    def toString1(self):
          stri = str(self.debutI) +"-"+str(self.milieu1I)+"-"+str(self.milieu2I) +"-"+str(self.finI)
          return stri

    def addSpace(self,n):
          s =""
          for i in range(0,n):
                  s += "  "
          return s

    def shape(self):
          stri = ""
          for h in self.helices:
                  stri += "("+ h.shape()+")"
          return stri
        
    def fillTheGap(self,n1,n2):
          tab = []
          for n  in range(n1,n2+1):
                  tab.append(n)
          return tab
        
    def createJ(self,jonction1):
          self.junction1 = jonction1
          debutI = self.milieu1I
          finI = self.milieu2I
          helices = self.helices
          #print ("len (self.helices) : "+str(len (self.helices)))
          if(len (self.helices) > 0 ):
                  allNt = self.fillTheGap(self.milieu1I,self.helices[0].debutI)
                  for i in range(0,len(self.helices)-1):
                          allNt = allNt + self.fillTheGap(self.helices[i].finI,self.helices[i+1].debutI)
                  allNt = allNt + self.fillTheGap(self.helices[-1].finI,self.milieu2I)
                  self.junction2 =  Junction(debutI,finI,helices,allNt,self.indice)
          else:
                  self.junction2 =   Junction(debutI,finI,helices,self.fillTheGap(debutI,finI),self.indice)
    

    def toStringForJS(self):
          stri = "{\"indice\":" + str(self.indice)+","
          stri += "\"pair\":["
          pairStrArray = []
          #print "self.pairTab : "+str(self.pairTab)
          for pair in self.pairTab:
                  pairStrArray.append(self.pairToString(pair))
          stri += ",".join(pairStrArray)
          stri += "],\"junction1\":"
          if self.junction1 is not None:
                  stri += self.jonction1.toStringForJS()
          else:
                  stri += "null"
          stri += ",\"junction2\":" + self.junction2.toStringForJS()
          
          stri += "}"
          return stri
        
    def pairToString(self,pair):
          stri = "["
          stri += str(pair[0])+","
          stri += str(pair[1])
          stri += "]"
          return stri
        
        
    
class Junction:
  count = 0
  def __init__(self,debutI,finI, helices,allNt,j):
          self.indice = j
          self.debutI = debutI
          self.finI = finI
          self.allNt = allNt
          self.helices = helices
          self.poigne = []
          for h in self.helices:
                  h.createJ(self)
                  self.poigne.append(h.pairTab[0])
          
  def printAllNt(self):
          print (";".join(str(i) for i in self.allNt))
    
  def toStringForJS(self):
          stri = "{\"indice\":" + str(self.indice) +" ,"
          stri += "\"poignes\":["
          poigneTab = []
          for po in self.poigne:
                  poigneTab.append(self.poigneToString(po))
          stri += ",".join(poigneTab)
          stri += "],\"allNt\":["
          #print "AllNt : "+",".join([str(x) for x in self.allNt])
          stri += ",".join([str(x) for x in self.allNt])
          stri += "],\"helices\":["
          heliceStrTab = []
          for h in self.helices:
                  heliceStrTab.append(h.toStringForJS())
          stri += ",".join(heliceStrTab)
          stri += "]}"
          return stri
    
  def poigneToString(self,pair):
          stri = "["
          stri += str(pair[0])+","
          stri += str(pair[1])
          stri += "]"
          return stri 

#iparens = iter('()')
#parens = dict(zip(iparens, iparens))
#closing = parens.values()

def balanced(astr):
          stack = []
          for c in astr:
                  d = parens.get(c, None)
                  if d:
                          stack.append(d)
                  elif c in closing:
                          if not stack or c != stack.pop():
                                  return False
          return not stack

def test():
  rna2d = RNA2Daio("gggggccccc")
  rna2d.displayRNA2D()

def window(fseq):
  w_size=5
  diAns = {}
  l = list(fseq)
  for i in range(0,len(l)-w_size+1):
          tempTab = []
          for i2 in range(i,i+w_size):
                  tempTab.append(l[i2])
          w = "".join(tempTab)
          if(w in diAns):
                  diAns[w]+=1
          else:
                  diAns[w]=1
  return diAns

def getMaskMcff(seq,di):
  di_ans={}
  w_size=5
  for w in di:
          st = list('x'*len(seq))
          indexTab = [m.start() for m in re.finditer(w,seq)]
          #print "indexTab"+str(indexTab)
          for i in indexTab:
                  for i2 in range(i,(i+w_size)):
                          st[i2] = '.'
          di_ans[w]="".join(st)
  return di_ans

def getMaskRso(seq,di,w_size=5):
  di_ans = {}
  for w in di:
          st = list('.'*len(seq))
          indexTab = [m.start() for m in re.finditer(w,seq)]
          for i in indexTab:
              for i2 in range(i,i+w_size):
                  st[i2] = 'x'
          di_ans[w]="".join(st)
  return di_ans
  
  
def fillString(s,ma):
  while(len(s) < ma):
          s = s + '-'
  return s