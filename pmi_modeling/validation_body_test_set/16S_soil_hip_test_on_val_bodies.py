#!/usr/bin/env python
# coding: utf-8

# ### Packages

# In[1]:


from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing, svm, metrics
from sklearn.model_selection import LeaveOneGroupOut, KFold, GridSearchCV, GroupKFold, cross_val_score, train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.datasets import make_regression
from scipy.stats import randint as sp_randint
from matplotlib import pyplot as plt
from collections import defaultdict


# In[2]:


import joblib
import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')
import biom
import calour as ca
from scipy import stats
import scipy
import pickle
import time
import math
import inspect
import operator
pd.set_option('display.max_rows', 100)


# In[3]:


def shared_features(test_dataset, training_dataset):
    def convert_to_df(dataset):
        matrix = dataset.data.todense()
        feature = list(dataset.feature_metadata.index)
        sample = list(dataset.sample_metadata.index)
        df = pd.DataFrame(matrix, columns=feature, index=sample)
        return df
    test_df = convert_to_df(test_dataset)
    training_df = convert_to_df(training_dataset)
    feature_ids = test_df.columns.intersection(training_df.columns)
    return feature_ids


# ### Import validation data and filter to keep soil samples

# In[4]:


validation = ca.read_amplicon('all-validation-normalized-species-table.biom','all-validation-bodies-metadata.txt', min_reads=1, normalize = None)
validation = validation.filter_samples('model_type', 'soil')
validation = validation.filter_prevalence(cutoff=1, fraction=0)
print(validation.sample_metadata.body_id.value_counts())
print(validation.sample_metadata.sample_type.value_counts())
print(validation.feature_metadata.count())


# ### Load original model data

# In[5]:


pmi = ca.read_amplicon('species-normalized-table.biom', 'combined-metadata-simple-oct2020.txt', min_reads=1, normalize = None)


# In[6]:


pmi = pmi.filter_samples('control', 'n')
pmi = pmi.filter_samples('soil_control', 'n')
pmi = pmi.filter_samples('sample_site', 'soil.hip')
pmi = pmi.filter_prevalence(cutoff=1, fraction=0)
print(pmi.sample_metadata.host_subject_id.value_counts())
print(pmi.sample_metadata.sample_site.value_counts())
print(pmi.feature_metadata.count())


# ### make validation and training sets match in features

# In[7]:


feature_ids = shared_features(pmi,validation)
pmi_match = pmi.filter_ids(feature_ids)
print(pmi_match.feature_metadata.count())
validation_match = validation.filter_ids(feature_ids)
print(validation_match.feature_metadata.count())


# ### build model

# In[8]:


X = pmi_match.data
y = pmi_match.sample_metadata['add_0c']
y = (y.astype(float))
regressor = RandomForestRegressor(n_estimators=1000, random_state=999, criterion='mae', max_depth=None, bootstrap=False, max_features=0.2, n_jobs=-1)
regressor.fit(X, y)
joblib.dump(regressor, 'PMI3_all_16S_soil_hip_species_general.pkl')
regressor = joblib.load('PMI3_all_16S_soil_hip_species_general.pkl')


# ### Testing model

# In[9]:


valX = validation_match.data
valy = validation_match.sample_metadata['add_0c']
valy = (valy.astype(float))


# In[10]:


yhat = regressor.predict(valX)


# In[11]:


mean_absolute_error(valy, yhat)


# In[12]:


results = pd.DataFrame(valy)
results['predicted_add'] = np.array(yhat)
results.index.name = None
f = open('soil_hip_validation_results.csv', 'w')
print("SampleID,true_add,predicted_add",file=f)
print(results.to_csv(header=False), file=f)


# In[13]:


results


# ## testing model with randomly assigned ADDs

# In[14]:


valy_rand = validation_match.sample_metadata['rand_add_0c']
valy_rand = (valy_rand.astype(float))
mean_absolute_error(valy_rand, yhat)


# In[15]:


results_rand = pd.DataFrame(valy_rand)
results_rand['predicted_add'] = np.array(yhat)
results_rand.index.name = None
f = open('soil_hip_validation_rand_results.csv', 'w')
print("SampleID,rand_add,predicted_add",file=f)
print(results_rand.to_csv(header=False), file=f)


# In[16]:


results_rand


# In[ ]:




