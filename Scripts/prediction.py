#########################################################################
#               PREDICT THE TOPOLOGY OF A GIVEN DATASET                 #
#########################################################################

import predictionmodules
slidingsize = 31


#Read the dataset to predict and save it in binary mode into a dictionary
dictionary2,filelines2 = predictionmodules.dicttopredict("testfile.txt")


#Get the input vector for topology prediction
testinginput = predictionmodules.getinputtesting(dictionary2,slidingsize)


#Predict the topology of the input dataset
predictedtopology = predictionmodules.predicting("finalizedmodel3.pkl",testinginput,dictionary2)


#Printing the right format for the predicted topologies (ID, seq, predicted topology)
predictionmodules.outputformat(filelines2,predictedtopology)




