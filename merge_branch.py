# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:45:52 2016

@author: Michaelhua
"""

import re
import csv
import operator


bpBranch={}
mfBranch={}
ccBranch={}


for astrLine in csv.reader( open("/Users/Yulei/Desktop/complete_GO_table.csv"), csv.excel_tab ):
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



#k=9
#right=float(0)
#total=float(0)
#for astrLine in csv.reader( open("D:/HSPH/BIO 508/Final/GONet/total_neighbours.csv"), csv.excel_tab ):    
#    bpGoList={}
#    #ccGoList={}
#    #mfGoList={}
#    geneName=astrLine[0]
#    #neighbors=[]
#    goLabels=[]
#    trueLabel=bpBranch.get(geneName)
#    if(trueLabel==None):
#        continue
#    else:
#        trueLabels=trueLabel.split(",")
#    
#    #wrong=float(0)
#
#    for j in range(1,k+1):
#        neighbor=astrLine[j]
#        if(bpBranch.get(neighbor)!=None):
#            goLabels=bpBranch.get(neighbor).split(",")
#            for g in goLabels:
#                if (bpGoList.get(g)!=None):
#                    bpGoList[g]+=1
#                else:
#                    bpGoList[g]=1
#            
#    predLabel=max(bpGoList, key=bpGoList.get)
#    total+=1
#    for label in trueLabels:
#        if(predLabel==label):
#            right+=1
#    outLine=geneName+"\t"+predLabel+"\t"+trueLabel+"\n"
#    outFile.write(outLine)
#
#print float(right/total)  
#outFile.close()

#outFile_bpCount=open("/Users/Yulei/Desktop/testCount_bp.txt",'w')
#outFile_ccCount=open("/Users/Yulei/Desktop/testCount_cc.txt",'w')
#outFile_mfCount=open("/Users/Yulei/Desktop/testCount_mf.txt",'w')

outFile_bp=open("/Users/Yulei/Desktop/accuracy_bp.txt",'w')
outFile_cc=open("/Users/Yulei/Desktop/accuracy_cc.txt",'w')
outFile_mf=open("/Users/Yulei/Desktop/accuracy_mf.txt",'w')
outFile_all=open("/Users/Yulei/Desktop/accuracy_all.txt",'w')

k=9

right_bp=float(0)
right_cc=float(0)
right_mf=float(0)

total_bp=float(0)
total_cc=float(0)
total_mf=float(0)

bpNone=float(0)
ccNone=float(0)
mfNone=float(0)

bpBranchCount=float(0)
ccBranchCount=float(0)
mfBranchCount=float(0)


for k in range(9, 104, 5):

    for astrLine in csv.reader( open("/Users/Yulei/Desktop/total_neighbours.csv"), csv.excel_tab ):    
        
        bpGoList={}
        ccGoList={}
        mfGoList={}

        geneName=astrLine[0]
        #neighbors=[]
        goLabels=[]

        trueLabel_bp=bpBranch.get(geneName)
        trueLabel_mf=mfBranch.get(geneName)
        trueLabel_cc=ccBranch.get(geneName)

        #print "trueLabel_bp: ", trueLabel_bp, 
        #print "trueLabel_cc: ", trueLabel_cc
        #print "trueLabel_mf: ", trueLabel_mf

        if(trueLabel_bp!=None):
            trueLabels_bp=trueLabel_bp.split(",") # i.e. trueLabels_bp: GO:0008152, GO:0009987

            bpBranchCount=float(0)
            ccBranchCount=float(0)
            mfBranchCount=float(0)

            for j in range(1,k+1):
                neighbor=astrLine[j]
                if(bpBranch.get(neighbor)!=None):
                    goLabels=bpBranch.get(neighbor).split(",") # get lables from neighbors, count frequency
                    bpBranchCount+=len(goLabels)
                    for g in goLabels:
                        if (bpGoList.get(g)!=None):
                            bpGoList[g]+=1
                        else:
                            bpGoList[g]=1

            if bpGoList!={}:
                predLabel_bp=max(bpGoList, key=bpGoList.get) # use neighbor with the largest count as prediction
                total_bp+=1

            for label in trueLabels_bp:
                if(predLabel_bp==label):
                    right_bp+=1

            #outFile_bpCount.write(str(bpBranchCount)+"\n")


        else: #if gene doesn't in the branch, check if neighbors don't give predicition of the branch, either
            bpNone+=1

            bpCount=float(0)
            ccCount=float(0)
            mfCount=float(0)

            for j in range(1,k+1):
                neighbor=astrLine[j]
                if(bpBranch.get(neighbor)!=None):
                    goLabels=bpBranch.get(neighbor).split(",")
                    bpCount+=len(goLabels)

        

              
        if(trueLabel_cc!=None):
            trueLabels_cc=trueLabel_cc.split(",")

            bpBranchCount=float(0)
            ccBranchCount=float(0)
            mfBranchCount=float(0)

            for j in range(1,k+1):
                neighbor=astrLine[j]
                if(ccBranch.get(neighbor)!=None):
                    goLabels=ccBranch.get(neighbor).split(",")
                    ccBranchCount+=len(goLabels)
                    for g in goLabels:
                        if (ccGoList.get(g)!=None):
                            ccGoList[g]+=1
                        else:
                            ccGoList[g]=1

            if ccGoList!={}:
                predLabel_cc=max(ccGoList, key=ccGoList.get)
                total_cc+=1

            for label in trueLabels_cc:
                if(predLabel_cc==label):
                    right_cc+=1

            #outFile_ccCount.write(str(ccBranchCount)+"\n")

        else: #if gene doesn't in the branch, check if neighbors don't give predicition of the branch, either
            ccNone+=1

            bpCount=float(0)
            ccCount=float(0)
            mfCount=float(0)

            for j in range(1,k+1):
                neighbor=astrLine[j]
                if(ccBranch.get(neighbor)!=None):
                    goLabels=ccBranch.get(neighbor).split(",")
                    ccCount+=len(goLabels)

            


        if(trueLabel_mf!=None):
            trueLabels_mf=trueLabel_mf.split(",")

            bpBranchCount=float(0)
            ccBranchCount=float(0)
            mfBranchCount=float(0)

            for j in range(1,k+1):
                neighbor=astrLine[j]
                if(mfBranch.get(neighbor)!=None):
                    goLabels=mfBranch.get(neighbor).split(",")
                    mfBranchCount+=len(goLabels)
                    for g in goLabels:
                        if (mfGoList.get(g)!=None):
                            mfGoList[g]+=1
                        else:
                            mfGoList[g]=1
            
            if mfGoList!={}:
                predLabel_mf=max(mfGoList, key=mfGoList.get)
                total_mf+=1

            for label in trueLabels_mf:       
                if(predLabel_mf==label):
                    right_mf+=1

            #outFile_mfCount.write(str(mfBranchCount)+"\n")


        else: #if gene doesn't in the branch, check if neighbors don't give predicition of the branch, either
            mfNone+=1

            bpCount=float(0)
            ccCount=float(0)
            mfCount=float(0)

            for j in range(1,k+1):
                neighbor=astrLine[j]
                if(mfBranch.get(neighbor)!=None):
                    goLabels=mfBranch.get(neighbor).split(",")
                    mfCount+=len(goLabels)



    #Accuracy for predicting the branch
    outLine_bp=str(k) + "\t" + str( float(right_bp/total_bp) )+ "\n"
    outLine_cc=str(k) + "\t" + str( float(right_cc/total_cc) )+ "\n"
    outLine_mf=str(k) + "\t" + str( float(right_mf/total_mf) )+ "\n"
    outLine_all=str(k) + "\t" + str( float( (right_bp+right_cc+right_mf)/(total_bp+total_cc+total_mf) ) ) + "\n"


    #Accuracy for predicting the GoList 
 
    
    outFile_bp.write(outLine_bp)
    outFile_cc.write(outLine_cc)
    outFile_mf.write(outLine_mf)
    outFile_all.write(outLine_all)



outFile_bp.close()    
outFile_cc.close()  
outFile_mf.close()    
outFile_all.close()  





