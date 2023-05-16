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


allShotgun = ca.read_amplicon('../MAG_species_table_final.biom','../shotgun_metadata_q2.txt', min_reads=1, normalize = None)

# ## Remove controls and samples below threshold
allShotgun = allShotgun.filter_samples('control', 'n')
allShotgun = allShotgun.filter_samples('below_threshold', 'n')
allShotgun.sample_metadata.description.value_counts()

# # Split by sampling location (soil v. skin)

# ## Soil sample filtering

Shotgun_soil = allShotgun.filter_samples('sample_type', 'soil')
Shotgun_soil = Shotgun_soil.filter_abundance(min_abundance=1)
Shotgun_soil.sample_metadata.sample_type.value_counts()

print("Number of samples: ",len(Shotgun_soil.sample_metadata.shotgun_barcode_used.value_counts()))

print("Number of bodies: ",len(Shotgun_soil.sample_metadata.host_subject_id.value_counts()))

# # Soil General Model using ADD 0C as Response and Adding Metadata as Features

Shotgun_soil_meta = Shotgun_soil.add_sample_metadata_as_features(['facility'])
print("Number of features: ",len(Shotgun_soil_meta.feature_metadata))

# designate data
X = Shotgun_soil_meta.data
y = Shotgun_soil_meta.sample_metadata['add_0c']
y = (y.astype(float))

# group by body
groups = Shotgun_soil_meta.sample_metadata['host_subject_id']

# fit regressor
final_regressor = RandomForestRegressor(n_estimators=1000, random_state=999, criterion='mae', bootstrap=False, max_depth=None, max_features=0.2, n_jobs=-1)
final_regressor.fit(X, y)
joblib.dump(final_regressor, '../models/PMI3_MAG_species_meta_general.pkl')
model_MAG_species_meta = joblib.load('../models/PMI3_MAG_species_meta_general.pkl')

# Determine important features
importances = model_MAG_species_meta.feature_importances_
std = np.std([tree.feature_importances_ for tree in model_MAG_species_meta.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
feature_metadata = Shotgun_soil_meta.feature_metadata

#print the important ids in order
general_importances = []
count = 0
print("Feature:\t\t\t\tImportance:")
for i in indices:
   general_importances += (Shotgun_soil_meta.feature_metadata.index.values[i], importances[indices[count]])
   if count < 5:
       print(str(count+1)+". "+str(Shotgun_soil_meta.feature_metadata.index.values[i])+"\t"+str(importances[indices[count]]))
   count += 1

print("Number of features present:", int(len(general_importances)/2))
general_importances_df = pd.DataFrame(np.array(general_importances).reshape(int(len(general_importances)/2),2))
np.savetxt("general_MAG_species_importances_meta.csv", general_importances_df, delimiter=",", fmt='%s')
