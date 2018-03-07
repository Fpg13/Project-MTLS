#For 9th of March deadline:

---- All the files mentioned below are contained in the Scripts folder ----

I have submitted the tweaked version of the predictor for this week. This time I have 4 important scripts:
-	One script for training the model: training.py
-	One script for the modules of the training script: trainingmodules.py
-	One script for predicting a topology: prediction.py
-	One script for the modules of the prediction script: predictionmodules.py

The parameters for SVM training are:
-	Sliding window size = 31
-	Kernel = Linear
-	C = 0.1

The trained model is saved as “finalizedmodel3.pkl”. This model has been trained with 16 sequences out of the 42 sequences of the whole dataset, since the size of the model trained on the whole dataset was 77MB.

You can test the topologies contained in “testfile.txt” with the script “prediction.py”

