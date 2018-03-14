#PREDICTION MODULES

import trainingmodules
import itertools
from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.externals import joblib
from sklearn import model_selection
import collections
import matplotlib.pyplot as plt



#################################################################
#       SAVE IN A BINARY DICTIONARY THE DATASET TO PREDICT      #
#################################################################

def dicttopredict(filetopredict):

    binaryaadict=trainingmodules.binarydict()
    dictionary2=collections.OrderedDict() 
    
    filehandle2=open(filetopredict,"r")
    filelines2=filehandle2.read().splitlines() 

    for i in range(len(filelines2)):
        if filelines2[i].startswith(">"):
            seq=[]
            seq=(filelines2[i+1])
            dictionary2[filelines2[i]]=seq
 
    for protein in dictionary2.keys():   
        binaryaalist2=[] 
    
        for residue in dictionary2[protein]:
            binaryaalist2.append(binaryaadict[residue])
        dictionary2[protein]=binaryaalist2 
    
    return(dictionary2,filelines2)
    
    filehandle2.close()
    
#################################################################
#               FUNCTION TO GET TESTING INPUT                   #
#################################################################
    
def getinputtesting(dictionary2,slidingsize):

    flankingaa=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    amountflanking=slidingsize//2
    aatoadd=[]
   
    for i in range(amountflanking):
        aatoadd=aatoadd+[flankingaa]
    
    modifieddictionary2=collections.OrderedDict() 

    for proteins in dictionary2.keys():
        modifieddictionary2[proteins]=aatoadd+dictionary2[proteins]+aatoadd

    testsetaa=[] #testing input
        
    for proteins in modifieddictionary2.keys():
        sequence=modifieddictionary2.get(proteins)
        for i in range(len(sequence)-2*amountflanking):
            tosave=sequence[i:(i+slidingsize)]
            merged2 = list(itertools.chain(*tosave))
            testsetaa.append(merged2)

    return(testsetaa)
    
    
#################################################################
#             LOAD THE MODEL AND PREDICT TOPOLOGY               #
#################################################################
        
def predicting(modelfilename,testsetaa,dictionary2):
    
    loaded_model=joblib.load(modelfilename)
    predictedtopology=loaded_model.predict(testsetaa)
 
    #Convert topology numbers back to letters:
    topologytype=[]
    for numbers in predictedtopology:
        if numbers==0:
            topologytype.append("i")
        elif numbers==1:
            topologytype.append("M")
        elif numbers==2:
            topologytype.append("o")
        
    return(topologytype)
    
#################################################################
#          OUTPUT TOPOLOGIES ALONG WITH ID AND SEQUENCE         #
#################################################################

def outputformat(filelines2,topologytype):

    initlen=0
    finallen=0
    for i in range(len(filelines2)):
        if filelines2[i].startswith(">"):
            print(filelines2[i])
            print(filelines2[i+1])
            finallen=finallen+len(filelines2[i+1])
            x ="".join(topologytype[initlen:finallen])
            print(x)
            initlen=finallen

 
#################################################################
#              CALCULATE F1 SCORE OF THE PREDICTION             #
################################################################# 
 
def f1score(topologytype,originalfile):
    filehandle=open(originalfile,"r")
    filelines=filehandle.read().splitlines() 
    truetopologies=[]
    topologytypelist=list(topologytype)
    for i in range(len(filelines)):
        if filelines[i].startswith(">"):
            truetopologies.extend(filelines[i+2])
    score=f1_score(truetopologies,topologytypelist,average="weighted")
    print("The F1 score for this prediction is:",score)
    
 
 
 
 
#################################################################
#                        CALL FUNCTIONS                         #
#################################################################

if __name__ == "__main__":
    
    print("")


