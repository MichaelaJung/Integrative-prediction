# -*- coding: utf-8 -*-

################################################################################
###                                                                          ###
### SHAP computation from Deep Learning models                               ###
###                                                                          ###
### by Carles Quesada-Traver, 2024                                           ###
###                                                                          ###
################################################################################
"""
This script uses the previously generated models to estimate the Shapely additive
explanations (SHAP) and stores them in a numpy format for further processing

GPU computation is advised, otherwise computation times might be excessive.

Usage: Run the script for each trait with the trait as the only argument. 
(i.e. "python F_SHAP_computation.py Acidity")
"""
#%%
import tensorflow as tf
import shap 
import numpy as np
import pandas as pd
import os
import sys
from sklearn.preprocessing import MinMaxScaler

#%%
#Load the data and define what goes in 
trait = sys.argv[1]
save_path = "Output/DeepLearning/Interpretation/"+ trait
if not os.path.exists(save_path):
    os.makedirs(save_path)
    
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
cv_template = pd.read_csv("Input/CV_Strategy/"+trait+"_CV.csv")
cv_template["ID_template"] = cv_template[['Envir','Genotype']].apply("-".join, axis = 1)
cv_template.set_index('ID_template', inplace = True)
pheno_subset_trait['ID_pheno'] = pheno_subset_trait[['Envir','Genotype']].apply("-".join, axis = 1)
pheno_subset_trait.set_index('ID_pheno', inplace = True)
cv_template = cv_template.reindex(pheno_subset_trait.index)
#Subset the response variable (trait)
y = pheno_subset_trait[trait].to_numpy()
y = y.reshape((y.shape[0],1)) 

#%% Scaling the data with the MinMaxScaler in the range (0,1) suitable for training

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


#%% Iterate over splits to compute and save SHAP values 

splits = cv_template.columns[7:]
for split in splits:
    snp_all = np.load(data_numpy_dir+"/snps/SNP_"+ trait + "_"+ split +".npy").astype('float32')
    #Split the different datasets between train and test
    test_index = cv_template[split] == 1
    train_index = cv_template[split] == 0
    weather_train,weather_test = weather_sc[train_index], weather_sc[test_index]
    pca_train, pca_test = pca_sc[train_index], pca_sc[test_index]
    snp_train, snp_test = snp_all[train_index], snp_all[test_index]
    soil_train, soil_test = soil_sc[train_index], soil_sc[test_index]
    y_train, y_test = y_sc[train_index],  y_sc[test_index]
    
    # Load the model
    model = tf.keras.models.load_model("Output/DeepLearning/Models/"+ trait +
                                       "/model_" + split +".keras")
    
    #Calculate SHAP values with Gradient Explainer and divide them for storage
    explainer = shap.GradientExplainer(model, [weather_train, pca_train, snp_train, soil_train])
    shap_values = explainer.shap_values([weather_test, pca_test, snp_test, soil_test])
    weather_SHAP = shap_values[0][0]
    pca_SHAP = shap_values[0][1]
    snp_SHAP = shap_values[0][2]
    soil_SHAP = shap_values[0][3]
    
    #Save them
    np.save(save_path +"/" +split+"_weather.npy",weather_SHAP)
    np.save(save_path +"/"+split+"_PCA.npy",pca_SHAP)
    np.save(save_path +"/"+split+"_SNP.npy",snp_SHAP)
    np.save(save_path +"/"+split+"_soil.npy",soil_SHAP)
