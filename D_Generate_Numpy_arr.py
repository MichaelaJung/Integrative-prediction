# -*- coding: utf-8 -*-
################################################################################
###                                                                          ###
### Generate Numpy arrays                                                    ###
###                                                                          ###
### by Carles Quesada-Traver, 2024                                           ###
###                                                                          ###
################################################################################

"""
In this script numpy arrays with the same number of row as the pheno file are
generated and saved. These will constitute the inmediate input for the deep 
learning models generated in the following script

Usage: Run the script (i.e. "python D_Generate_Numpy_arr.py")
"""

# Script to generate ready to use data formats for training the DL architecture
import pandas as pd
import numpy as np

# Directory where numpy data will be saved
data_dir = "Output/Intermediate/numpy_arrays"

# Add important columns to the data for further weather processing
data = pd.read_csv("Output/Intermediate/Pheno_processed_pivoted.csv")
data['Year'] = data.Envir.str.split(".").str[1]
data['END'] = data['Year'] + '-11-1'
data['END'] = pd.to_datetime(data['END'], infer_datetime_format = True)
data['START'] = data['Year'] + '-01-01'
data['START'] = pd.to_datetime(data['START'], infer_datetime_format = True)
#%% Load weather data and transform dates for easier operability
all_weather = pd.read_csv("Output/Intermediate/Weather_daily.csv")
all_weather['Envir'] = all_weather['Location']+'.' + all_weather['Date'].str.split("-").str[0]
all_weather['Date'] = pd.to_datetime(all_weather['Date'], infer_datetime_format= True)


# Define variables controlling the dimensions of the output 3D object
n_days = (data.iloc[0]['END'] - data.iloc[0]['START']).days - 4 # 300 days
len_interval =  1 # Change if longer intervals want to be averaged for input
n_samples = data.shape[0]
n_variables = 3   
n_timesteps = int(n_days/len_interval) 

#Create an empty numpy array that will be the output of the processing 
X_weather = np.zeros((n_samples, n_timesteps, n_variables)) 

for x in range(n_samples):
    # Filter for the evironment of that sample
    env_mask = pd.array(all_weather['Envir'] == data['Envir'][x])
    env_df = all_weather[env_mask]
    # Its better to fill the dataframe from 1st of January 
    mask_dates = pd.array(env_df['Date'] < data['END'][x])
    # Obtain the following n_days from 1st of January
    df = env_df[mask_dates].tail(n_days)
    df_np = df.to_numpy()
    data_reshaped = np.zeros((n_timesteps, n_variables)) 
    for i in range(0,n_timesteps):
	    range_1 = i * len_interval
        range_2 = (i + 1) * len_interval
        data_to_avg = df_np[range_1:range_2, :]        
        data_reshaped[i, 0] = np.mean(data_to_avg[:,2],axis=0) # Temperature / Avged
        data_reshaped[i, 1] = np.mean(data_to_avg[:,3],axis=0) # Humidity / Avged
        data_reshaped[i, 2] = np.sum(data_to_avg[:,4],axis=0) # Radiation / summed, change to average if needed
    
    X_weather[x,:,:] = data_reshaped

np.save("%s/weather.npy"%(data_dir), X_weather)
#%% PCA common for every model // General genetic variation 

PCA = pd.read_csv("Output/Intermediate/PCAir_PCs.csv")
PCA = PCA.rename(columns={ PCA.columns[0]: "Genotype" })
n_comp = PCA.shape[1] - 1
X_PCA = np.zeros((n_samples,n_comp))

for x in range(n_samples):
    row = PCA[PCA['Genotype'] == data['Genotype'][x]].to_numpy()[:,1:]
    X_PCA[x,:] = row
    
np.save("%s/pca.npy"%(data_dir), X_PCA)
#%% Phenotypes to filter missing data per trait and save SNP data 

pheno = data.drop(['Year', 'END', 'START'], axis=1)
traits = pheno.columns[2:]
for trait in traits:
    splits = pd.read_csv("Output/Intermediate/CV_Strategy/"+trait+"_CV.csv").columns[2:]
    # Filter the phenotype rows based on missing values 
    nan_mask_trait = np.isnan(pheno[trait])
    pheno_subset_trait = pheno[~nan_mask_trait].reset_index(drop=True)
    n_samples_snp = pheno_subset_trait.shape[0]
    for split in splits: 
        snp_trait_split = pd.read_csv("Output/Intermediate/geno_files/"+trait+"/geno_"+split+".csv")
        snp_trait_split = snp_trait_split.rename(columns={ snp_trait_split.columns[0]: "Genotype" })
        # We remove the "G_" Added to make it suitable for modeling before
        snp_trait_split['Genotype'] = snp_trait_split['Genotype'].str.slice(2)
        # Empy numpy array for n_samples for each trait
        n_snps = snp_trait_split.shape[1]-1
        X_SNP_trait = np.zeros((n_samples_snp, n_snps))
        
        for x in range(n_samples_snp):
            row = snp_trait_split[snp_trait_split['Genotype'] == pheno_subset_trait['Genotype'][x]].to_numpy()[:,1:]
            X_SNP_trait[x,:] = row
        np.save( data_dir + '/snps/SNP_'+ trait + '_'+split+'.npy', X_SNP_trait)

#%%## Soil 
soil_data = pd.read_csv("Output/Intermediate/Soil.csv")
n_soil_parameters = soil_data.shape[1] - 1
X_Soil = np.zeros((n_samples,n_soil_parameters))

for x in range(n_samples): 
    row = soil_data[soil_data['Envir'] == data['Envir'][x]].to_numpy()[:,1:]
    X_Soil[x,:] = row
    
np.save("%s/soil.npy"%(data_dir), X_Soil)





