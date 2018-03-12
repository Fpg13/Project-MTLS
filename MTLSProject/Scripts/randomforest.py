#########################################################################
#                              RANDOM FOREST                            #
#########################################################################


import trainingmodules
import probando
from sklearn.ensemble import RandomForestRegressor


# Instantiate model with 1000 decision trees
rf=RandomForestRegressor(n_estimators = 1000, random_state = 42)

svminput,svmoutput=probando.pssmconvert()

clf=rf.fit(svminput, svmoutput)
modelfilename="rfmodel.pkl"
joblib.dump(clf,modelfilename)
print("The model was successfully saved!")
