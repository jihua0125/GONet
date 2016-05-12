# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:15:06 2016

@author: Michaelhua
"""

import re
import csv
import operator


bpBranch={}
mfBranch={}
ccBranch={}


for astrLine in csv.reader( open("D:/HSPH/BIO 508/Final/GONet/complete_GO_table.csv"), csv.excel_tab ):
    goType=astrLine[2]  
    geneName=astrLine[0]
    goLabel=astrLine[1]
    if (goType=="bp"):
        if(bpBranch.get(geneName)==None):
            bpBranch[geneName]=goLabel
        else:
            bpBranch[geneName]=bpBranch[geneName]+","+goLabel
            
    elif (goType=="mf"):
        if(mfBranch.get(geneName)==None):
            mfBranch[geneName]=goLabel
        else:
            mfBranch[geneName]=mfBranch[geneName]+","+goLabel
    
    else:
        if(ccBranch.get(geneName)==None):
            ccBranch[geneName]=goLabel
        else:
            ccBranch[geneName]=ccBranch[geneName]+","+goLabel
            
predFile=open("D:/HSPH/BIO 508/Final/GONet/prediction.csv",'w')           
k=9

for astrLine in csv.reader( open("D:/HSPH/BIO 508/Final/GONet/unknown_total_neighbours.csv"), csv.excel_tab ):    
    
    bpGoList={}
    ccGoList={}
    mfGoList={}

    geneName=astrLine[0]
    goLabels=[]


    for j in range(1,k+1):
        neighbor=astrLine[j]
        if(bpBranch.get(neighbor)!=None):
            goLabels=bpBranch.get(neighbor).split(",") # get lables from neighbors, count frequency
            for g in goLabels:
                if (bpGoList.get(g)!=None):
                    bpGoList[g]+=1
                else:
                    bpGoList[g]=1

    if bpGoList!={}:
        predLabel_bp=max(bpGoList, key=bpGoList.get) # use neighbor with the largest count as prediction
    else:
        predLabel="NA"
    outline=geneName+"\t"+predLabel_bp+"\tbp\n"
    predFile.write(outline)

    #outFile_bpCount.write(str(bpBranchCount)+"\n")    

          
    for j in range(1,k+1):
        neighbor=astrLine[j]
        if(ccBranch.get(neighbor)!=None):
            goLabels=ccBranch.get(neighbor).split(",") # get lables from neighbors, count frequency
            for g in goLabels:
                if (ccGoList.get(g)!=None):
                    ccGoList[g]+=1
                else:
                    ccGoList[g]=1

    if ccGoList!={}:
        predLabel_cc=max(ccGoList, key=ccGoList.get) # use neighbor with the largest count as prediction
    else:
        predLabel="NA"
    outline=geneName+"\t"+predLabel_cc+"\tcc\n"
    predFile.write(outline)
    
    for j in range(1,k+1):
        neighbor=astrLine[j]
        if(mfBranch.get(neighbor)!=None):
            goLabels=mfBranch.get(neighbor).split(",") # get lables from neighbors, count frequency
            for g in goLabels:
                if (mfGoList.get(g)!=None):
                    mfGoList[g]+=1
                else:
                    mfGoList[g]=1

    if mfGoList!={}:
        predLabel_mf=max(mfGoList, key=mfGoList.get) # use neighbor with the largest count as prediction
    else:
        predLabel="NA"
    outline=geneName+"\t"+predLabel_mf+"\tmf\n"
    predFile.write(outline)
    predFile.write("\n")

predFile.close()    
    


