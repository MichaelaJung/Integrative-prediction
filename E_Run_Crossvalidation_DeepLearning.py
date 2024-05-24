# -*- coding: utf-8 -*-

################################################################################
###                                                                          ###
###  Genomic prediction following a crossvalidation scheme, using the deep   ###
###  learning model                                                          ###
###                                                                          ###
###   by Carles Quesada-Traver, 2024                                         ###
###                                                                          ###
################################################################################

"""
This script runs the crossvalidation within the Deep Learning framework for 
a given trait. It takes as input the previously generated numpy arrays, and
produces 3 types of outputs, the models, the predictions and the general metrics.

GPU computation is advised, otherwise computation times might be excessive.

Usage: Run the script for each trait with the trait as the only argument. 
(i.e. "python E_Run_Crossvalidation_DeepLearning.py Acidity")
"""

# Importing required packages 
# General packages
import os
import tensorflow as tf
from keras import backend as K
import numpy as np
import pandas as pd
import sys

# Important packages for model training
from deep_model import *
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import pearsonr
from math import sqrt
from keras.optimizers import Adam
import time

# Setting up GPU environment for compatibility between Tensorflow and Keras
os.environ["CUDA_VISIBLE_DEVICES"] = "0"  
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "1"
os.environ["KERAS_BACKEND"] = "tensorflow"
config = tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
sess = tf.compat.v1.Session()
K.set_session(sess)

#%% Define some important variables for training 
trait = sys.argv[1] # Run the program with trait name as argument
batch_size = 256 
epochs = 1500 # Enough epochs, that will be stopped according to training data 
lr = 0.0001 

#%% Load datasets shared among traits
data_numpy_dir = "Output/Intermediate/numpy_arrays"
weather_all = np.load(data_numpy_dir + "/weather.npy").astype('float32')
pca_all = np.load(data_numpy_dir + "/pca.npy").astype('float32')
soil_all = np.load(data_numpy_dir + "/soil.npy").astype('float32')

# Load the phenotypes and select the one of interest
pheno = pd.read_csv("Output/Intermediate/Pheno_processed_pivoted.csv")
# Filter for missing values
nan_mask_trait = np.isnan(pheno[trait])
pheno_subset_trait = pheno[~nan_mask_trait].reset_index(drop=True)

# Change CV template order to match phenotype order
cv_template = pd.read_csv("Output/Intermediate/CV_Strategy/"+trait+"_CV.csv")
cv_template["ID_template"] = cv_template[['Envir','Genotype']].apply("-".join, axis = 1)
cv_template.set_index('ID_template', inplace = True)
pheno_subset_trait['ID_pheno'] = pheno_subset_trait[['Envir','Genotype']].apply("-".join, axis = 1)
pheno_subset_trait.set_index('ID_pheno', inplace = True)
cv_template = cv_template.reindex(pheno_subset_trait.index)
#Subset the response variable (trait)
y = pheno_subset_trait[trait].to_numpy()
y = y.reshape((y.shape[0],1)) 

#%% Scaling the data with the MinMaxScaler in the range (0,1) suitable for training
from sklearn.preprocessing import MinMaxScaler

# For weather scaling we have to transform to a 2D, scale and come back to 3D
weather_trait = weather_all[~nan_mask_trait]
weather_2D = weather_trait.reshape(weather_trait.shape[0], weather_trait.shape[1] * weather_trait.shape[2])
scaler_weather = MinMaxScaler(feature_range = (0,1))
scaler_weather = scaler_weather.fit(weather_2D)
weather_2D_sc = scaler_weather.transform(weather_2D)
weather_sc = weather_2D_sc.reshape(weather_trait.shape[0], weather_trait.shape[1],weather_trait.shape[2])

## PCA, soil and response variables can be directly scaled
pca_trait = pca_all[~nan_mask_trait]
scaler_pca = MinMaxScaler(feature_range = (0,1))
scaler_pca = scaler_pca.fit(pca_trait)
pca_sc = scaler_pca.transform(pca_trait)

soil_trait = soil_all[~nan_mask_trait]
scaler_soil = MinMaxScaler(feature_range = (0,1))
scaler_soil  = scaler_soil.fit(soil_trait)
soil_sc = scaler_soil.transform(soil_trait)


scaler_y = MinMaxScaler(feature_range = (0,1))
scaler_y = scaler_y.fit(y)
y_sc = scaler_y.transform(y)

#%% Set up the variables to store the predictions
columns = ['Split','r', 'time', 'r2', 'RMSE', 'MAE']
df = pd.DataFrame(columns = columns)
save_name = "Split_snps"
save_dir = "Output/DeepLearning/"
# Create a copy of the CV_template to save all the predictions.
predictions_full = cv_template 
splits = cv_template.columns[2:]
# Early stopping callback for stopping training when training loss stops improving
stop_early = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=35, restore_best_weights= True)

#%% Iteratively perform crossvalidation
for split in splits:
    # Load the SNP data which is fold specific
    snp_all = np.load(data_numpy_dir+"/snps/SNP_"+ trait + "_"+ split +".npy").astype('float32')
    # Start computing the time 
    time_0 = time.time()
    tf.keras.backend.clear_session() # Clears previously defined variables
    print("Processing split %s" % split)

    
    ## Separate input data into train and test splits
    test_index = cv_template[split] == 1
    train_index = cv_template[split] == 0    
    weather_train,weather_test = weather_sc[train_index], weather_sc[test_index]
    pca_train, pca_test = pca_sc[train_index], pca_sc[test_index]
    snp_train, snp_test = snp_all[train_index], snp_all[test_index]
    soil_train, soil_test = soil_sc[train_index], soil_sc[test_index]
    y_train, y_test = y_sc[train_index],  y_sc[test_index]
    
    ## Initialize the model with the necessary variables
    time_steps = weather_train.shape[1]
    n_pca = pca_train.shape[1]
    n_snp = snp_train.shape[1]
    n_soil = soil_train.shape[1]
    pred_model = model(time_steps, n_pca, n_snp, n_soil)
    pred_model.compile(loss = 'mean_squared_error', optimizer = Adam(learning_rate = lr))
    # Initialize fitting
    hist = pred_model.fit([weather_train, pca_train, snp_train, soil_train], y_train,
                          batch_size = batch_size,
                          epochs = epochs,
                          verbose = 2, 
                          shuffle = False,
                          callbacks = [stop_early])
    # Perform prediction of test instances
    y_pred = pred_model.predict([weather_test, pca_test, snp_test, soil_test])
    # Reshape and bring back to initial values 
    y_pred = scaler_y.inverse_transform(y_pred).reshape(y_test.shape[0],)
    # Record time needed for prediction
    time_1 = time.time()
    elapsed_time = time_1 - time_0
    # Save the model for further analysis
    modelname = "model_" + str(split)+ ".keras"
    save_model_path = save_dir +"Models/"+ trait +"/"+ modelname
    os.makedirs(save_model_path, exist_ok=True) 
    pred_model.save(save_model_path)
    y_true = scaler_y.inverse_transform(y_test).reshape(y_test.shape[0],)
    
    # Compute interesting metrics
    rmse = sqrt(mean_squared_error(y_true, y_pred))
    mae = sqrt(mean_absolute_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    r = pearsonr(y_true, y_pred)[0]
    row = ["Split %s" % split, r, elapsed_time, r2, rmse, mae]
    df.loc[len(df)] = row
    
    # Compute full predictions and store in respective dataframe
    y_pred_all = pred_model.predict([weather_sc, pca_sc, snp_all, soil_sc])
    y_pred_all = scaler_y.inverse_transform(y_pred_all).reshape(y_pred_all.shape[0],)
    predictions_full[split] = y_pred_all 

# Compute and save general metrics    
mean_rmse = df['RMSE'].mean() 
mean_mae = df['MAE'].mean()
mean_r = df['r'].mean()
mean_t = df['time'].mean()
mean_r2 = df['r2'].mean()
means = ["Mean", mean_r, mean_t, mean_r2, mean_rmse, mean_mae]
df.loc[len(df)] = means
df.to_csv(save_dir + trait +"_metrics_splits.csv", index = False)


#%% Format the results to match the other modelling strategy output

# Melt the CV_template 
df_melted_cv = cv_template.melt(id_vars=['Envir', 'Genotype'], 
                    var_name='CV', 
                    value_name='Testing')
# Extract the run and fold information
df_melted_cv['Run'] = df_melted_cv['CV'].str.extract(r'CV(\d)_Split')[0]
df_melted_cv['Fold'] = df_melted_cv['CV'].str.extract(r'Split(\d)')[0]

# Convert Run and Fold to numeric types
df_melted_cv['Run'] = pd.to_numeric(df_melted_cv['Run'])
df_melted_cv['Fold'] = pd.to_numeric(df_melted_cv['Fold'])

# Add a column of observed values to the predictions
df_observed = pheno_subset_trait[['Envir', 'Genotype', trait]]
df_observed.columns = ['Envir', 'Genotype', 'Observed']
predictions_full_obs = pd.merge(predictions_full, df_observed, how = 'left')

# Melt the predictions with the observed values to match previous data
df_melted_pred_obs = predictions_full_obs.melt(id_vars=['Envir', 'Genotype', 'Observed'], 
                    var_name='CV', 
                    value_name='Predicted')

# Append the columns of observed and predicted phenotypic values
df_melted_cv['Predicted'] = df_melted_pred_obs['Predicted']
df_melted_cv['Observed'] = df_melted_pred_obs['Observed']


# Drop the CV column as it is no longer needed and save
df_result = df_melted_cv.drop(columns=['CV'])
df_result.to_csv(save_dir+"Predictions/"+trait+"_predictions_formated_csv", index = False)