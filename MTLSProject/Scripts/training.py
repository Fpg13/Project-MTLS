#########################################################################
#                 TRAINING OF A MODEL FROM INPUT DATASET               #
#########################################################################

import trainingmodules


#Open the file and store it into a dictionary
#dictionary = trainingmodules.saveindict("example.txt")
dictionary = trainingmodules.saveindict("membrane-beta_3state.3line.txt")


#Convert the dictionary into binary values and get a list with all the topologies (training output for svm)
binarydictionary, trainingoutput = trainingmodules.convertbinary(dictionary)


#Add flanking vectors, compute sliding windows and merge them into a single list for svm input
slidingsize = 31
traininginput = trainingmodules.getinputsvm(dictionary,slidingsize)


#Train the model and get the accuracy. Then save it
scores, clf, savedmodel = trainingmodules.svmtraining(traininginput,trainingoutput)
print(scores)
print("The predictor accuracy with a sliding window of %s is: %0.5f (+/+ %0.5f)" % (slidingsize,        scores.mean(), scores.std()*2)) 



