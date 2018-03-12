#TRAINING MODULES

import itertools
from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
from sklearn.externals import joblib
from sklearn import model_selection
import collections
import matplotlib.pyplot as plt


#################################################################
#   FUNCTION TO SAVE ID, SEQUENCE AND TOPOLOGY IN A DICTIONARY  #
#################################################################

def saveindict(file):

    dictionary=collections.OrderedDict() #The keys are the protein IDs, and the values are a list [aa seq, topology]
    filehandle=open(file,"r")
    filelines=filehandle.read().splitlines() 
    
    for i in range(len(filelines)):
        if filelines[i].startswith(">"):
            #filelines[i]=filelines[i].replace(">","")
            seqandtopology=[0,0]
            seqandtopology[0]=(filelines[i+1])
            seqandtopology[1]=(filelines[i+2])
            dictionary[filelines[i]]=seqandtopology

    return(dictionary)

    filehandle.close()
    
#################################################################
#            FUNCTION TO CREATE BINARY DICTIONARY               #
#################################################################

def binarydict():

    binaryaadict={} #The keys are the amino acid sequences and the values are the binary codes.

    binaryaadict["A"]=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["C"]=[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["D"]=[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["E"]=[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["F"]=[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["G"]=[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["H"]=[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["I"]=[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["K"]=[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["L"]=[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0]
    binaryaadict["M"]=[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]
    binaryaadict["N"]=[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
    binaryaadict["P"]=[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]
    binaryaadict["Q"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]
    binaryaadict["R"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0]
    binaryaadict["S"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0]
    binaryaadict["T"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]
    binaryaadict["V"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0]
    binaryaadict["W"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0]
    binaryaadict["Y"]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]   

    return(binaryaadict)


#################################################################
#            FUNCTION TO CONVERT DATASET INTO BINARY            #
#################################################################

def convertbinary(dictionary):

    converteddictionary=dictionary
    trainingsettopology1=[] #set for training output of svm.
    binaryaadict=binarydict()
    
    for protein in converteddictionary.keys():   
        binaryaalist=[] #List that stores the aa sequence in binary code.
        topologybinary=[] #List that stores the topology in binary code.
    
        for residue in converteddictionary[protein][0]:
            binaryaalist.append(binaryaadict[residue])
        converteddictionary[protein][0]=binaryaalist 
    
        for letter in converteddictionary[protein][1]:
            if letter=="i":
                topologyi=0
                topologybinary.append(topologyi)
            elif letter=="M":
                topologyM=1
                topologybinary.append(topologyM)
            elif letter=="o":
                topologyo=2
                topologybinary.append(topologyo)
        converteddictionary[protein][1]=topologybinary 
    
        trainingsettopology1.append(topologybinary)

    trainingsettopology=[]
    for sublist in trainingsettopology1:
        for item in sublist:
            trainingsettopology.append(item)

    return(converteddictionary,trainingsettopology)


#################################################################
#            FUNCTION TO GET TRAINING INPUT TO SVM              #
#################################################################

def getinputsvm(dictionary,slidingsize):

    flankingaa=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    amountflanking=slidingsize//2
    aatoadd=[]

    for i in range(amountflanking):
        aatoadd=aatoadd+[flankingaa]
        modifieddictionary=collections.OrderedDict() #flanking aa also stored, but without the topology.

    for proteins in dictionary.keys():
        modifieddictionary[proteins]=aatoadd+dictionary[proteins][0]+aatoadd

    #Storing sliding windows into lists

    trainingsetaa=[] #training input

    for proteins in modifieddictionary.keys():
        sequence=modifieddictionary.get(proteins)
        for i in range(len(sequence)-2*amountflanking):
            tosave=sequence[i:(i+slidingsize)]
            merged = list(itertools.chain(*tosave))
            trainingsetaa.append(merged)

    return(trainingsetaa)
        
    
#################################################################
#       CROSS VALIDATED SETS & SVM TRAINING & SAVE MODEL        #
#################################################################   

def svmtraining(trainingsetaa,trainingsettopology):
    
    clf=svm.SVC(kernel="linear",C=0.1,gamma=0.01,cache_size=3000).fit(trainingsetaa,trainingsettopology)
    scores=cross_val_score(clf,trainingsetaa,trainingsettopology,cv=5,scoring="accuracy",n_jobs=-1) 


    ##CONFUSION MATRIX##
    
    predicted=cross_val_predict(clf,trainingsetaa,trainingsettopology,cv=5)
    conf_mat=confusion_matrix(trainingsettopology,predicted)
    print(conf_mat)
    plt.plot(conf_mat)
    plt.title("Confusion matrix")
    plt.xlabel("Predicted label")
    plt.ylabel("True label")
    plt.show()
   

    modelfilename="PSSMmodel.pkl"
    joblib.dump(clf,modelfilename)
    print("The model was successfully saved!")
    
    return(scores,clf,modelfilename)


#################################################################
#                        CALL FUNCTIONS                         #
#################################################################

if __name__ == "__main__":
    
    print("")


    
