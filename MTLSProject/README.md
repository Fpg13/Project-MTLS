--------------------------
FINAL PREDICTOR SUBMISSION
--------------------------

In a few words:

The main script you have to use is “prediction.py”. From here you can choose different FILES TO PREDICT, as well as load different SAVED MODELS and calculate the F1 score of the prediction by selecting the FILE that contains the REAL TOPOLOGIES:

- “testfile.txt”: the proteins in this file can be predicted with the linear model (“finalizedmodel.pkl”), the random forest model (“rfmodel.pkl”) and the decision tree model (“decisiontreemodel.pkl”). To get the F1 score of this prediction, use “completetestfile.txt”, which contains the real topologies of these proteins.
- “fastadataset.fasta”: the proteins in this file (50 new extracted proteins) can be predicted with the linear model (“finalizedmodel.pkl”), the random forest model (“rfmodel.pkl”) and the decision tree model (“decisiontreemodel.pkl”). To get the F1 score of this prediction, use “modifiednewdataset.txt”, which contains the real topologies of these 50 new extracted proteins.
- “testPSSM.txt”:  the PSSM files of the proteins in this file have been previously obtained by running the psi-blast bash script and their topologies can be predicted with the PSSM model (“PSSMmodel.pkl”). To get the F1 score of this prediction, use “completefilePSSM.txt”, which contains the real topologies of these proteins.


Extra clarification of the scripts/models/text files:

- "training.py" calls the functions contained in "trainingmodules.py", necessary to train each of the models.
- "prediction.py" calls the functions contained in "predictionmodules.py", necessary for topology prediction of different files with each of the selected models.
- There are 3 linear models: "finalizedmodel.pkl" has been created with training of the whole dataset "membrane-beta_3state.3line.txt"; "finalizedmodel2.pkl" has been created with training of the 21 sequences of the original dataset, contained in "example.txt"; "finalizedmodel3.pkl" has been created with training of the 16 sequences of the original dataset, contained in "example2.txt". This was done because the size of the first 2 models was too big.
- "PSSMmodel.pkl" is created with the script in "pssmconverter.py".
- "prepareinputpssm.py" prepares the input data so that it is compatible with the PSSM model.
- "rfmodel.pkl" is created with the script in "randomforest.py".
- "decisiontreemodel.pkl" is created with the script in "decisiontree.py".
- "psiblast.py" is used to create one FASTA file for each of the proteins in the original dataset.
- The psiblast and pssm files contained in the Datasets folder have been obtained with the bash script "psiblast.sh".
- "newdataset.txt" contains the 50 extra proteins extracted from TOPDB."modifynewdataset.py" adjusts the topologies of these proteins so that it is compatible with the script.

