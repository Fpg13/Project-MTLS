#########################################################################
#              PREPARE INPUT FOR PREDICTION WITH PSSM MODEL             #
#########################################################################

import collections
import os
import itertools

def convertPSSMinput(PSSMfile):
    svminput=[]
    slidingsize=31
    dictionary=collections.OrderedDict() #The keys are the protein IDs, and the values are a list [aa seq, topology]
    filehandle=open(PSSMfile,"r")
    filelines=filehandle.read().splitlines() 
    
    for i in range(len(filelines)):
        if filelines[i].startswith(">"):
            seq=[]
            seq.append(filelines[i+1])
            dictionary[filelines[i]]=seq

    for proteins in dictionary.keys():
        if os.path.isfile("../Datasets/PSSMfiles/"+proteins+".fasta.pssm"):		
            print("Running:"+proteins+".fasta.pssm")		
            filehandle=open("../Datasets/PSSMfiles/"+proteins+".fasta.pssm", 'r')
            filelines=filehandle.read().splitlines()
            flankingaa=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
            amountflanking=slidingsize//2
            aatoadd=[]
            for i in range(amountflanking):
                aatoadd=[flankingaa]+aatoadd
            for i in range(3,len(filelines)-6):
                scores=filelines[i].split()
                aatoadd=aatoadd+[(list(scores[22:-2]))]
            for i in range(amountflanking):
                aatoadd=aatoadd+[flankingaa]
        
            #Normalization step#
            for locations in aatoadd:
                for number in range(0,len(locations)):
                    locations[number]=int(locations[number])/100
            dictionary[proteins][0]=aatoadd

            for i in range(len(aatoadd)-2*amountflanking):
                tosave=aatoadd[i:(i+slidingsize)]
                merged=list(itertools.chain(*tosave))
                svminput.append(merged)

    return(svminput)

