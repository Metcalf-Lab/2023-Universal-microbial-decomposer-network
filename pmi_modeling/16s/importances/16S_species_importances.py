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


###########################################
####### SKIN FACE
###########################################
# Importing data, no TSS normalization performed here since table is already normalized
all_16S = ca.read_amplicon('../bioms/species-normalized-table.biom','../combined-metadata-simple-oct2020.txt', min_reads=1, normalize = None)

# Remove controls
all_16S = all_16S.filter_samples('control', 'n')
all_16S = all_16S.filter_samples('soil_control', 'n')

# sample filtering
all_16S_skin_face = all_16S.filter_samples('sample_site', 'skin.face')
all_16S_skin_face = all_16S_skin_face.filter_abundance(min_abundance=1)
all_16S_skin_face_meta = all_16S_skin_face.add_sample_metadata_as_features(['facility'])

# designate data
X = all_16S_skin_face_meta.data
y = all_16S_skin_face_meta.sample_metadata['add_0c']
y = (y.astype(float))

# group by body
groups = all_16S_skin_face_meta.sample_metadata['host_subject_id']

# fix the sparsity issue
type(X)
X.data = np.nan_to_num(X.data)

# fit regressor
final_regressor = RandomForestRegressor(n_estimators=1000, random_state=999, criterion='mae', bootstrap=False, max_depth=None, max_features=0.2, n_jobs=-1)
final_regressor.fit(X, y)
joblib.dump(final_regressor, '../models/PMI3_all_16S_skin_face_meta_general.pkl')
model_all_16S_skin_face_meta = joblib.load('../models/PMI3_all_16S_skin_face_meta_general.pkl')

# Determine important features
importances = model_all_16S_skin_face_meta.feature_importances_
std = np.std([tree.feature_importances_ for tree in model_all_16S_skin_face_meta.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
feature_metadata = all_16S_skin_face_meta.feature_metadata

#print the important ids in order
general_importances = []
count = 0
print("Feature:\t\t\t\tImportance:")
for i in indices:
   general_importances += (all_16S_skin_face_meta.feature_metadata.index.values[i], importances[indices[count]])
   if count < 5:
       print(str(count+1)+". "+str(all_16S_skin_face_meta.feature_metadata.index.values[i])+"\t"+str(importances[indices[count]]))
   count += 1

print("Number of features present:", int(len(general_importances)/2))
general_importances_df = pd.DataFrame(np.array(general_importances).reshape(int(len(general_importances)/2),2))
np.savetxt("general_16S_skin_face_importances_meta.csv", general_importances_df, delimiter=",", fmt='%s')

###########################################
####### SKIN HIP
###########################################
# Importing data, no TSS normalization performed here since table is already normalized
all_16S = ca.read_amplicon('../bioms/species-normalized-table.biom','../combined-metadata-simple-oct2020.txt', min_reads=1, normalize = None)

# Remove controls
all_16S = all_16S.filter_samples('control', 'n')
all_16S = all_16S.filter_samples('soil_control', 'n')

# sample filtering
all_16S_skin_hip = all_16S.filter_samples('sample_site', 'skin.hip')
all_16S_skin_hip = all_16S_skin_hip.filter_abundance(min_abundance=1)
all_16S_skin_hip_meta = all_16S_skin_hip.add_sample_metadata_as_features(['facility'])

# designate data
X = all_16S_skin_hip_meta.data
y = all_16S_skin_hip_meta.sample_metadata['add_0c']
y = (y.astype(float))

# group by body
groups = all_16S_skin_hip_meta.sample_metadata['host_subject_id']

# fix the sparsity issue
type(X)
X.data = np.nan_to_num(X.data)

# fit regressor
final_regressor = RandomForestRegressor(n_estimators=1000, random_state=999, criterion='mae', bootstrap=False, max_depth=None, max_features=0.2, n_jobs=-1)
final_regressor.fit(X, y)
joblib.dump(final_regressor, '../models/PMI3_all_16S_skin_hip_meta_general.pkl')
model_all_16S_skin_hip_meta = joblib.load('../models/PMI3_all_16S_skin_hip_meta_general.pkl')

# Determine important features
importances = model_all_16S_skin_hip_meta.feature_importances_
std = np.std([tree.feature_importances_ for tree in model_all_16S_skin_hip_meta.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
feature_metadata = all_16S_skin_hip_meta.feature_metadata

#print the important ids in order
general_importances = []
count = 0
print("Feature:\t\t\t\tImportance:")
for i in indices:
   general_importances += (all_16S_skin_hip_meta.feature_metadata.index.values[i], importances[indices[count]])
   if count < 5:
       print(str(count+1)+". "+str(all_16S_skin_hip_meta.feature_metadata.index.values[i])+"\t"+str(importances[indices[count]]))
   count += 1

print("Number of features present:", int(len(general_importances)/2))
general_importances_df = pd.DataFrame(np.array(general_importances).reshape(int(len(general_importances)/2),2))
np.savetxt("general_16S_skin_hip_importances_meta.csv", general_importances_df, delimiter=",", fmt='%s')

###########################################
####### SOIL HIP
###########################################
# Importing data, no TSS normalization performed here since table is already normalized
all_16S = ca.read_amplicon('../bioms/species-normalized-table.biom','../combined-metadata-simple-oct2020.txt', min_reads=1, normalize = None)

# Remove controls
all_16S = all_16S.filter_samples('control', 'n')
all_16S = all_16S.filter_samples('soil_control', 'n')

# sample filtering
all_16S_soil_hip = all_16S.filter_samples('sample_site', 'soil.hip')
all_16S_soil_hip = all_16S_soil_hip.filter_abundance(min_abundance=1)
all_16S_soil_hip_meta = all_16S_soil_hip.add_sample_metadata_as_features(['facility'])

# designate data
X = all_16S_soil_hip_meta.data
y = all_16S_soil_hip_meta.sample_metadata['add_0c']
y = (y.astype(float))

# group by body
groups = all_16S_soil_hip_meta.sample_metadata['host_subject_id']

# fix the sparsity issue
type(X)
X.data = np.nan_to_num(X.data)

# fit regressor
final_regressor = RandomForestRegressor(n_estimators=1000, random_state=999, criterion='mae', bootstrap=False, max_depth=None, max_features=0.2, n_jobs=-1)
final_regressor.fit(X, y)
joblib.dump(final_regressor, '../models/PMI3_all_16S_soil_hip_meta_general.pkl')
model_all_16S_soil_hip_meta = joblib.load('../models/PMI3_all_16S_soil_hip_meta_general.pkl')

# Determine important features
importances = model_all_16S_soil_hip_meta.feature_importances_
std = np.std([tree.feature_importances_ for tree in model_all_16S_soil_hip_meta.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
feature_metadata = all_16S_soil_hip_meta.feature_metadata

#print the important ids in order
general_importances = []
count = 0
print("Feature:\t\t\t\tImportance:")
for i in indices:
   general_importances += (all_16S_soil_hip_meta.feature_metadata.index.values[i], importances[indices[count]])
   if count < 5:
       print(str(count+1)+". "+str(all_16S_soil_hip_meta.feature_metadata.index.values[i])+"\t"+str(importances[indices[count]]))
   count += 1

print("Number of features present:", int(len(general_importances)/2))
general_importances_df = pd.DataFrame(np.array(general_importances).reshape(int(len(general_importances)/2),2))
np.savetxt("general_16S_soil_hip_importances_meta.csv", general_importances_df, delimiter=",", fmt='%s')

###########################################
####### SOIL FACE
###########################################
# Importing data, no TSS normalization performed here since table is already normalized
all_16S = ca.read_amplicon('../bioms/species-normalized-table.biom','../combined-metadata-simple-oct2020.txt', min_reads=1, normalize = None)

# Remove controls
all_16S = all_16S.filter_samples('control', 'n')
all_16S = all_16S.filter_samples('soil_control', 'n')

# sample filtering
all_16S_soil_face = all_16S.filter_samples('sample_site', 'soil.face')
all_16S_soil_face = all_16S_soil_face.filter_abundance(min_abundance=1)
all_16S_soil_face_meta = all_16S_soil_face.add_sample_metadata_as_features(['facility'])

# designate data
X = all_16S_soil_face_meta.data
y = all_16S_soil_face_meta.sample_metadata['add_0c']
y = (y.astype(float))

# group by body
groups = all_16S_soil_face_meta.sample_metadata['host_subject_id']

# fix the sparsity issue
type(X)
X.data = np.nan_to_num(X.data)

# fit regressor
final_regressor = RandomForestRegressor(n_estimators=1000, random_state=999, criterion='mae', bootstrap=False, max_depth=None, max_features=0.2, n_jobs=-1)
final_regressor.fit(X, y)
joblib.dump(final_regressor, '../models/PMI3_all_16S_soil_face_meta_general.pkl')
model_all_16S_soil_face_meta = joblib.load('../models/PMI3_all_16S_soil_face_meta_general.pkl')

# Determine important features
importances = model_all_16S_soil_face_meta.feature_importances_
std = np.std([tree.feature_importances_ for tree in model_all_16S_soil_face_meta.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
feature_metadata = all_16S_soil_face_meta.feature_metadata

#print the important ids in order
general_importances = []
count = 0
print("Feature:\t\t\t\tImportance:")
for i in indices:
   general_importances += (all_16S_soil_face_meta.feature_metadata.index.values[i], importances[indices[count]])
   if count < 5:
       print(str(count+1)+". "+str(all_16S_soil_face_meta.feature_metadata.index.values[i])+"\t"+str(importances[indices[count]]))
   count += 1

print("Number of features present:", int(len(general_importances)/2))
general_importances_df = pd.DataFrame(np.array(general_importances).reshape(int(len(general_importances)/2),2))
np.savetxt("general_16S_soil_face_importances_meta.csv", general_importances_df, delimiter=",", fmt='%s')
