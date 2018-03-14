#########################################################################
#                            DECISION TREE                              #
#########################################################################


from datetime import datetime
from sklearn.externals import joblib
import trainingmodules
import pssmconverter 
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix

svminput,svmoutput=pssmconverter.pssmconvert()
s=3

tree=DecisionTreeClassifier(min_samples_split=s)
clf=tree.fit(svminput,svmoutput)
    
scores=cross_val_score(clf,svminput,svmoutput,cv=5,scoring="accuracy",n_jobs=-1) 
            
print("Scores for min_samples_split =",s,"are:",scores)
print("The decision tree predictor accuracy is: %0.5f (+/+ %0.5f)" % (scores.mean(), scores.std()*2)) 
    
modelfilename="decisiontreemodel.pkl"
joblib.dump(clf,modelfilename)
print("The model was successfully saved!")

#Confusion matrix DT model
predicted=cross_val_predict(clf,svminput,svmoutput,cv=5)
conf_mat=confusion_matrix(svmoutput,predicted)
print(conf_mat)
