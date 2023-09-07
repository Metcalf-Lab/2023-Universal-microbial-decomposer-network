#!/usr/bin/env python
# coding: utf-8

from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing, svm, metrics
from sklearn.model_selection import LeaveOneGroupOut, KFold, GridSearchCV, GroupKFold, cross_val_score, train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
#from sklearn.externals import joblib
from sklearn.datasets import make_regression
from scipy.stats import randint as sp_randint
from collections import defaultdict

import joblib
import pandas as pd
import numpy as np
import biom
import calour as ca
from scipy import stats
import scipy
import pickle
import time
import math
import inspect
import operator
pd.set_option('display.max_rows', 10000)


#############################################################################################################
# Load the metadata file
metadata_file = 'combined-metadata.txt'
metadata = pd.read_csv(metadata_file, sep='\t')

# Extract the 'megyesi_score' as X and ADD as y
X = metadata['megyesi_score']
unique_values_X = X.unique()
#print(unique_values_X)
X = (X.astype(int))

y = metadata['add_0c']
unique_values_y = y.unique()
#print(unique_values_y)
y = y.values + 1e-6
y = (y.astype(float))

print("X shape:", X.shape)
print("y shape:", y.shape)

# group by body
groups = metadata['host_subject_id']

# outer_cv creates 36 folds by leave-one-body-out for estimating generalization error
outer_cv = LeaveOneGroupOut().split(X, y, groups=groups)
# 
# prints the number of folds in the outer loop
print("Number of outer folds to perform: ", LeaveOneGroupOut().get_n_splits(X, y, groups=groups))
# 
# hyperparameter grid to test
param_grid = {"max_depth": [None, 4],
              "max_features": ['auto', 0.2],
              "bootstrap": [True, False]}
# 
# creates loop iteration counter and empty lists for storage
count=1
model_parameters = []
nested_cv_scores = []
# 
# loops through the sets of training and test ids in the outer loop
# the number of loops will match the number of folds as each fold is used
with open('TBS_allpred_output.csv', 'w') as f:
	print("SampleID,true_add,predicted_add",file=f)
	for train_ids, test_ids in outer_cv:
#     inner_cv creates 35 folds by leave-one-body-out for hyperparamter tuning
#     uses only the train ids present in the current outer loop fold which is one less body since the outer loop
#     folds are also using leave-one-body-out
		inner_cv = LeaveOneGroupOut().split(X[train_ids], y[train_ids], groups=groups[train_ids])
		print("X train shape:", X[train_ids].shape)
		print("y train shape:", y[train_ids].shape)
		print("group train shape:", groups[train_ids].shape)
#     setting rf parameters
		rf = RandomForestRegressor(n_estimators=500, random_state=999, criterion='mae')
		print("rf")
#     grid search cv using rf and the hyperparameter grid on the inner_cv training set
		rf_grid = GridSearchCV(estimator=rf, param_grid=param_grid, cv=inner_cv, n_jobs=-1, scoring='neg_mean_absolute_error')
		print("rf_grid")
#     fit the grid search model on the inner_cv training set, which will tell us the best
#     parameters chosen by that inner cv
		rf_grid.fit(X[train_ids].values.reshape(-1, 1), y[train_ids])
		print("rf_gridfit")
#     converts best params dict to string to save it
		res = ",".join(("{}={}".format(*i) for i in rf_grid.best_params_.items()))
#     attaches each loops best params
		model_parameters.append(res)
#     prints outer fold number the loop is currently on
		print("Outer fold:",count)
#     prints number of inner folds in the outer loop (should be the same each time at 35)
		print("Number of inner folds:",LeaveOneGroupOut().get_n_splits(X[train_ids], y[train_ids], groups=groups[train_ids]))
#     prints best param and CV score from this inner loop set 
		print("Best params:",rf_grid.best_params_)
		print("Best CV score (MAE):",-rf_grid.best_score_)
#     uses the best model created from the inner loop to predict the outer loop body left out
		yhat = rf_grid.predict(X[test_ids].values.reshape(-1, 1))
		results = pd.DataFrame(y[test_ids])
		results['predicted_add'] = np.array(yhat)
		results.index.name = None
		print(results.to_csv(header=False), file=f)
		MAE = mean_absolute_error(y[test_ids], yhat)
		print("Prediction score (MAE):",MAE)
		nested_cv_scores.append(MAE)
		print("**********************")
		count+=1
f.close()
# prints the mean of the nested cv score (generalization error) of the models
print("\TBS Nested CV score (generalization error): " + str(np.mean(nested_cv_scores)))
#
 
# define function to merge 2 lists as paired tuples    
def merge(list1, list2): 
    merged_list = tuple(zip(list1, list2))  
    return merged_list 

# merge parameters to nested scores as tuples
merged_list = merge(model_parameters, nested_cv_scores)

# puts the paired tuple list into a dictionary
params_and_scores = defaultdict(list)
for k, v in merged_list:
    params_and_scores[k].append(v)


# prints dict of params and scores, should see each set of parameters chosen and their scores
print("All parameters and scores")
print(params_and_scores)

 
# create a dict to assign the mean score to each param set
avgDict = {}
for k,v in params_and_scores.items():
    avgDict[k] = sum(v)/ float(len(v))
print("All parameters with avg scores")
print(avgDict)


# prints the param set with the best mean error
print("Best params:",min(avgDict, key=avgDict.get))
print("Best mean score:",min(avgDict.values()))
