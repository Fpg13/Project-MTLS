#########################################################################
#                            RANDOM FOREST                              #
#########################################################################

from sklearn.externals import joblib
import trainingmodules
import pssmconverter 
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix


svminput,svmoutput=pssmconverter.pssmconvert()
n=330
s=3

rf=RandomForestClassifier(n_estimators = n, min_samples_split = s)
clf=rf.fit(svminput, svmoutput)
scores=cross_val_score(clf,svminput,svmoutput,cv=5,scoring="accuracy",n_jobs=-1) 
        
print("Scores for n_estimators =",n,"and min_samples_split =",s,"are:",scores)
print("The RF predictor accuracy is: %0.5f (+/+ %0.5f)" % (scores.mean(), scores.std()*2)) 

modelfilename="rfmodel.pkl"
joblib.dump(clf,modelfilename,compress=9)
print("The model was successfully saved!")
    
#Confusion matrix RF model
predicted=cross_val_predict(clf,svminput,svmoutput,cv=5)
conf_mat=confusion_matrix(svmoutput,predicted)
print(conf_mat)
     

