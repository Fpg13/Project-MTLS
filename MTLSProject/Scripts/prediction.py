#########################################################################
#               PREDICT THE TOPOLOGY OF A GIVEN DATASET                 #
#########################################################################

import predictionmodules
import prepareinputpssm
slidingsize = 31

#------------------------------------------------------------------------------------------------#

#Read the dataset to predict (it only has ID and seq) and save it in binary mode into a dictionary

#filetotest="testfile.txt"
#filetotest="fastadataset.fasta" #50 new extracted proteins
filetotest="testPSSM.txt"

dictionary2,filelines2 = predictionmodules.dicttopredict(filetotest)

#------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------#

#Get the input vector for topology prediction

testinginput = predictionmodules.getinputtesting(dictionary2,slidingsize)
if filetotest=="testPSSM.txt":
    testinginput=prepareinputpssm.convertPSSMinput("testPSSM.txt")
    
#------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------#    

#Predict the topology of the input dataset using different models:

#modeltoload="finalizedmodel.pkl"
modeltoload="PSSMmodel.pkl"
#modeltoload="rfmodel.pkl"
#modeltoload="decisiontreemodel.pkl"
 
predictedtopology = predictionmodules.predicting(modeltoload,testinginput,dictionary2)

#------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------#

#Printing the right format for the predicted topologies (ID, seq, predicted topology)

predictionmodules.outputformat(filelines2,predictedtopology)

#------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------#

#Calculate the F1 score of the prediction. Original file contains ID, seq and real topology.

#originalfile="completetestfile.txt" #Contains real topologies of "testfile.txt"
#originalfile="modifiednewdataset.txt" #Contains real topologies of "fastadataset.fasta"
originalfile="completefilePSSM.txt" #Contains real topologies of "testPSSM.txt"

predictionmodules.f1score(predictedtopology,originalfile)

#------------------------------------------------------------------------------------------------#

