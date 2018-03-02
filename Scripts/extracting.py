
        
        ###########################
        ### CREATE DICTIONARIES ###
        ###########################

dictionary={} #dictionary where the keys are the protein ids, and the values are a list [aa seq, topology]
#filehandle=open("membrane-beta_3state.3line.txt","r")
filehandle=open("example.txt","r")
filelines=filehandle.read().splitlines() 
    
#print(filelines)

for i in range(len(filelines)):
    if filelines[i].startswith(">"):
        filelines[i]=filelines[i].replace(">","")
        seqandtopology=[0,0]
        seqandtopology[0]=(filelines[i+1])
        seqandtopology[1]=(filelines[i+2])
        dictionary[filelines[i]]=seqandtopology
           
#print(dictionary)

binaryaadict={} #dictionary where the keys are the amino acids and the values are the binary codes.

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

#print(binaryaadict)    
    

trainingsettopology1=[] #set for training output of svm.

for protein in dictionary.keys():   
    binaryaalist=[] #List that stores the aa sequence in binary code.
    topologybinary=[] #List that stores the topology in binary code.
    
    for residue in dictionary[protein][0]:
        binaryaalist.append(binaryaadict[residue])
    dictionary[protein][0]=binaryaalist #Turn the aa sequences of the dictionary into binary and replace them.
    
    for letter in dictionary[protein][1]:
        if letter=="i":
            topologyi=0
            topologybinary.append(topologyi)
        elif letter=="M":
            topologyM=1
            topologybinary.append(topologyM)
        elif letter=="o":
            topologyo=2
            topologybinary.append(topologyo)
    dictionary[protein][1]=topologybinary #Turn the topologies of the dictionary into binary and replace them.
    
    trainingsettopology1.append(topologybinary)

#print(dictionary)

#print(trainingsettopology1)
trainingsettopology=[]
for sublist in trainingsettopology1:
    for item in sublist:
        trainingsettopology.append(item)

#print(trainingsettopology)
#print(len(trainingsettopology))
#print(len(dictionary))
#print(dictionary)


        #######################
        ### SLIDING WINDOWS ###
        #######################
        
flankingaa=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
  
for size in range(3,4,2): #so it only takes the 3 so far
    slidingsize=size
    amountflanking=size//2
    aatoadd=[]

    for i in range(amountflanking):
        aatoadd=aatoadd+[flankingaa]
    
    modifieddictionary={} #dictionary with the flanking aa also stored, but without the topology.

    for proteins in dictionary.keys():
        modifieddictionary[proteins]=aatoadd+dictionary[proteins][0]+aatoadd

#print(modifieddictionary)


        ##########################################
        ### STORING SLIDING WINDOWS INTO LISTS ###
        ##########################################
        
    trainingsetaa=[] #training input

    import itertools #necessary to merge the lists of aa contained in a sliding window into a single list.

    for proteins in modifieddictionary.keys():
        sequence=modifieddictionary.get(proteins)
        for i in range(len(sequence)-2*amountflanking):
            tosave=sequence[i:(i+slidingsize)]
            merged = list(itertools.chain(*tosave))
            trainingsetaa.append(merged)


    #print(trainingsetaa)
    #print(trainingsettopology)
    #print(len(trainingsetaa))
    #print(len(trainingsettopology))
        
       
        #############################################
        ### CROSS VALIDATED SETS & SVM TRAINING ###      
        #############################################
        
    from sklearn import svm
    from sklearn.model_selection import cross_val_score
    
    clf=svm.SVC(kernel="linear",C=1).fit(trainingsetaa,trainingsettopology)
    scores=cross_val_score(clf,trainingsetaa,trainingsettopology,cv=5) #"estimate accuracy of linear kernel svm on our input dataset by splitting the data, fitting a model and scoring 5 consecutive times (with different splits each time)"
    
    print(scores)
    print("The predictor accuracy with a sliding window of %s is: %0.5f (+/+ %0.5f)" % (slidingsize, scores.mean(), scores.std()*2)) 
    
    #testsetaa=predicttopology(filetopredict)
    #clf.predict(testsetaa)







        ##########################################################################################
        ### SCRIPT TO CONVERT THE INPUT FILE INTO THE NECESSARY FORMAT FOR TOPOLOGY PREDICTION ###
        ##########################################################################################
        
    #Calling the function that does this does not work. Ask and call the function instead of pasting everything here.
    
    
        ###########################
        ### CREATE DICTIONARIES ###
        ###########################

    dictionary2={} 
    filetopredict="testfile.txt"
    filehandle2=open(filetopredict,"r")
    filelines2=filehandle2.read().splitlines() 
    
    #print(filelines2)

    for i in range(len(filelines2)):
        if filelines2[i].startswith(">"):
            filelines2[i]=filelines2[i].replace(">","")
            seq=[]
            seq=(filelines2[i+1])
            dictionary2[filelines2[i]]=seq
           
    #print(dictionary2)
 
    for protein in dictionary2.keys():   
        binaryaalist2=[] 
    
        for residue in dictionary2[protein]:
            binaryaalist2.append(binaryaadict[residue])
        dictionary2[protein]=binaryaalist2 
    
    #print(dictionary2)

    

        #######################
        ### SLIDING WINDOWS ###
        #######################
        
        flankingaa=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
  
    for size in range(3,4,2): #so it only takes the 3 so far.
        slidingsize=size
        amountflanking=size//2
        aatoadd=[]

        for i in range(amountflanking):
            aatoadd=aatoadd+[flankingaa]
    
        modifieddictionary2={} 

        for proteins in dictionary2.keys():
            modifieddictionary2[proteins]=aatoadd+dictionary2[proteins]+aatoadd


    #print(modifieddictionary2)


    
        ##########################################
        ### STORING SLIDING WINDOWS INTO LISTS ###
        ##########################################
        
        testsetaa=[] #training input

        import itertools 
        
        for proteins in modifieddictionary2.keys():
            sequence=modifieddictionary2.get(proteins)
            #print(proteins)
            #print(sequence)
            for i in range(len(sequence)-2*amountflanking):
                tosave=sequence[i:(i+slidingsize)]
                merged2 = list(itertools.chain(*tosave))
                testsetaa.append(merged2)

        
        #print(testsetaa)
        #print(len(testsetaa))
        
        predictedtopology=clf.predict(testsetaa)
        print("The predicted topology is:", predictedtopology)
        print(len(predictedtopology))

        












        ######################
        ### SAVE THE MODEL ###     
        ######################
        
    #from sklearn.externals import joblib
    
    #modelfilename="finalizedmodel.pkl"
    #joblib.dump(clf,modelfilename)  #I need to store all the clf, not only the last one.


        ######################
        ### LOAD THE MODEL ###
        ######################
    
    #from sklearn import model_selection

    #inputsequence=open("testfile.txt","r")
    #loaded_model=joblib.load(modelfilename)
    #result=loaded_model.predict(inputsequence)
    #print(result)


filehandle.close()
filehandle2.close()
    


#Call functions

#if __name__ == "__main__":
    
    #predicttopology(filetopredict)
    
