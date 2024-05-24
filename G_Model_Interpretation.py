# -*- coding: utf-8 -*-

################################################################################
###                                                                          ###
### Model Interpretation through SHAP values                                 ###
###                                                                          ###
### by Carles Quesada-Traver, 2024                                           ###
###                                                                          ###
################################################################################

"""
This script uses the previously generated SHAP values and summarizes them in two 
different ways:

    (I) Stream contribution: Perc. contribution per split of each model stream
    (II) Mean SHAP contribution: Absolute mean SHAP values averaged across folds

Usage: Run the script (i.e. "python G_Model_Interpretation.py")
"""

import numpy as np
import pandas as pd
import os
import glob


# Empty df to save
stream_contrib_df = pd.DataFrame(columns = ["Trait", "Split", "PC_stream", "SNP_stream", "Soil_stream", "Weather_stream"])
mean_SHAP_contribution = pd.DataFrame(columns = ["Trait","Feature", "meanSHAP", "n_folds"])
# Load the phenotypes and select the one of interest
pheno = pd.read_csv("Output/Intermediate/Pheno_processed_pivoted.csv")
traits = pheno.columns[2:]
for trait in traits:
    nan_mask_trait = np.isnan(pheno[trait])
    pheno_subset_trait = pheno[~nan_mask_trait].reset_index(drop=True)
    
    # Change CV template order to match Phenotype and all the other data 
    cv_template = pd.read_csv("Output/Intermediate/CV_Strategy/"+trait+"_CV.csv")
    cv_template["ID_template"] = cv_template[['Envir','Genotype']].apply("-".join, axis = 1)
    cv_template.set_index('ID_template', inplace = True)
    pheno_subset_trait['ID_pheno'] = pheno_subset_trait[['Envir','Genotype']].apply("-".join, axis = 1)
    pheno_subset_trait.set_index('ID_pheno', inplace = True)
    cv_template = cv_template.reindex(pheno_subset_trait.index)
    #Create a master df per trait to store Split, Feature and SHAP
    columns_master_df = ["Split", "Feature", "SHAP"]
    master_df = pd.DataFrame(columns = columns_master_df)
    # Iterate over splits for a given trait
    splits = cv_template.columns[2:]
    data_dir = "Output/DeepLearning/Interpretation/" + trait
    for split in splits:
        split_df = pd.DataFrame(columns = columns_master_df)
        # For each split list the 4 SHAP files created previously
        pattern = split+"_*"
        files = glob.glob(os.path.join(data_dir, pattern))
        for file in files:
            # Extract the type of file the loop is dealing with
            dtype = file.split(split)[1][1:-4]
            temp_df = pd.DataFrame(columns = columns_master_df)
            
            # Depending on the file type the operations will be different
            if dtype == "weather":
                array = np.load(file, allow_pickle=True)
                
                # Process 300 days of temperature
                temp_SHAPs = np.mean(np.abs(array[:,:,0]),0)
                temp_labels = [dtype + "_temp_day_"+ str(n+1) for n in range(300)]
                
                # Process 300 days of humidity
                humidity_SHAPs = np.mean(np.abs(array[:,:,1]),0)
                humidity_labels = [dtype + "_hum_day_"+ str(n+1) for n in range(300)]
                
                # Process 300 days of radiation
                radiation_SHAPs = np.mean(np.abs(array[:,:,2]),0)
                radiation_labels = [dtype + "_rad_day_"+ str(n+1) for n in range(300)]
                # Create the absolute values
                abs_mean_array = np.concatenate((temp_SHAPs, humidity_SHAPs, radiation_SHAPs))
                temp_df['Feature'] = temp_labels + humidity_labels + radiation_labels
            
            elif dtype == "SNP":
                array = np.load(file, allow_pickle=True)
                abs_mean_array = np.mean(np.abs(array),0)
                # Annotate with respective SNP names
                temp_df['Feature'] = pd.read_csv("Output/Intermediate/geno_files"+ trait+"/geno_"+split+".csv").columns[1:]       
    
            elif dtype == "soil":
                array = np.load(file, allow_pickle=True)
                abs_mean_array = np.mean(np.abs(array),0)
                n_feat = abs_mean_array.shape[0]
                # Annotate with specific soil features
                temp_df['Feature'] = ["soil_"+ feature for feature in pd.read_csv("Output/Intermediate/Soil.csv").columns[1:]]       
            
            else:
                array = np.load(file)
                abs_mean_array = np.mean(np.abs(array),0)
                # Annotate with PC number
                temp_df['Feature'] = [dtype + "_"+ str(n+1) for n in range(abs_mean_array.shape[0])]
            
            temp_df['SHAP'] = abs_mean_array
            temp_df['Split'] = split
            # Concatente the temporary dataframe with the split dataframe 
            split_df = pd.concat([split_df, temp_df], ignore_index=True, axis=0)
        
        # Per split calculate relative contribution of each stream
        SHAP_sum = np.sum(split_df['SHAP'])
        contrib_PCA  = split_df[split_df['Feature'].str.startswith("PCA")]['SHAP'].sum()/SHAP_sum*100
        contrib_SNP = split_df[split_df['Feature'].str.startswith("AX")]['SHAP'].sum()/SHAP_sum*100
        contrib_weather = split_df[split_df['Feature'].str.startswith("weather")]['SHAP'].sum()/SHAP_sum*100
        contrib_soil = split_df[split_df['Feature'].str.startswith("soil")]['SHAP'].sum()/SHAP_sum*100
        row_contrib = [trait, split, contrib_PCA, contrib_SNP, contrib_soil, contrib_weather]
        # Annotate the stream dataframe with a row per split of each trait
        stream_contrib_df.loc[len(stream_contrib_df)] = row_contrib
        # Append split results of absolute mean SHAP values to previous splits
        master_df = pd.concat([master_df, split_df], ignore_index=True, axis=0)
    
    # Summarize the values per trait
    summary_df = master_df.groupby('Feature').agg(
        meanSHAP=('SHAP', 'mean'),
        sdSHAP=('SHAP', 'std'),
        n_folds=('SHAP', 'size')  # Interesting for SNPs not it every fold
    ).reset_index()
    
    summary_df['Trait'] = trait
    summary_df = summary_df.fillna(0)
    # Add the trait summarized outputs to previous summarized outputs
    mean_SHAP_contribution = pd.concat([mean_SHAP_contribution, summary_df], ignore_index = True, axis = 0)
    
# Save the dataframes 
stream_contrib_df.to_csv("Output/DeepLearning/Interpretation/stream_contribution.csv", index = False)
mean_SHAP_contribution.to_csv("Output/DeepLearning/Interpretation/mean_SHAP_contribution.csv", index = False)

