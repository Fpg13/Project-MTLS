#########################################################################
#                  PSSM DATA CONVERSION FOR MODEL TRAINING              #
#########################################################################

import itertools
import os
import trainingmodules
import collections

def pssmconvert():

    slidingsize=31
    dictionary = trainingmodules.saveindict("membrane-beta_3state.3line.txt")
    svminput=[]
    svmoutput=[]

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

            topologybinary=[]
            for letter in dictionary[proteins][1]:
                if letter=="i":
                    topologyi=0
                    topologybinary.append(topologyi)
                elif letter=="M":
                    topologyM=1
                    topologybinary.append(topologyM)
                elif letter=="o":
                    topologyo=2
                    topologybinary.append(topologyo)
            dictionary[proteins][1]=topologybinary 
            svmoutput.extend(topologybinary)


#########################################################################
#                      TRAIN MODEL WITH PSSM DATA                       #
#########################################################################    
    
    scores, clf, savedmodel = trainingmodules.svmtraining(svminput,svmoutput)
    print(scores)
    print("The PSSM predictor accuracy with a sliding window of %s is: %0.5f (+/+ %0.5f)" % (slidingsize,        scores.mean(), scores.std()*2))
    
    return(svminput,svmoutput)


if __name__ == "__main__":
    
    svminput,svmoutput = pssmconvert()
