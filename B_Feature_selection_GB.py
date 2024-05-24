# -*- coding: utf-8 -*-

################################################################################
###                                                                          ###
### Feature Selection using GB                                               ###
###                                                                          ###
### by Carles Quesada-Traver, 2024                                           ###
###                                                                          ###
################################################################################

"""
Using the output from "A_PreProcess_Data_1.R", this script performs a systematic
feature selection of SNPs for each combination of trait and split, it produces
a common output file contaning relevan information for the pipeline.

Usage: Run the script (i.e. "python B_Feature_selection_GB.py")
"""

# Load the required packages
import numpy as np
import pandas as pd
import os
from sklearn.ensemble import GradientBoostingRegressor

# Set the working directory
GB_working_dir ="Output/Intermediate/GB_feature_selection/"

# Declare an empty df where the results will be saved
importance_df = pd.DataFrame(columns = ["Trait","Split", "SNP", "importance"])

# Read geno file which is common for all traits and splits
geno = pd.read_csv(GB_working_dir+ "all.geno")

# Load traits and iterate over them
traits = pd.read_csv("Output/Intermediate/Pheno_processed_pivoted.csv").columns[2:]
for trait in traits:
    # Obtain the split files containing the phenotypic data
    files_splits = os.listdir(GB_working_dir + trait)
    
    for file in files_splits:
        # Load phenotypic data for a certain split (Only training instances)
        pheno_split = pd.read_csv(GB_working_dir + trait +"/"+file)
        
        # Filter the genotypic matrix for the genotypes of interest
        mask_geno = geno[geno.columns[0]].isin(pheno_split["Genotype"])
        geno_split = geno[mask_geno].reset_index(drop=True)
        
        # Set genotype to index 
        pheno_split.set_index('Genotype', inplace = True)
        geno_split.set_index(geno_split.columns[0], inplace = True)
        
        # Reindex pheno to make sure the instances match
        pheno_split = pheno_split.reindex(geno_split.index)
        
        # Transform data to numpy for model fitting
        pheno_np = np.ravel(pheno_split.to_numpy())
        geno_np = geno_split.to_numpy()
        
        # Declare and fit the model with 200 estimators
        print("Fitting Gradient Boosting regressor for ", file)
        reg = GradientBoostingRegressor(random_state = 0, n_estimators = 200)
        reg.fit(geno_np, pheno_np)
        
        # Extract feature importance and filter for important snps
        feature_importance = reg.feature_importances_
        onlyimportantlogic = feature_importance >  0
        n_important = np.sum(onlyimportantlogic)
        sorted_idx = np.argsort(feature_importance)[::-1]
        important_idx = sorted_idx[0:n_important]
        important_snps = geno_split.columns[important_idx].tolist()
        importance = feature_importance[important_idx].tolist()
        
        # Organize in a df the selected SNPs per split
        data_splits = {'Trait':trait, 'Split': file[:-4],
                       'SNP': important_snps, 'importance': importance} 
        split_df = pd.DataFrame(data_splits)
        
        # Append them to the previously declared dataframe
        importance_df = pd.concat([importance_df, split_df], ignore_index=True, axis=0)
        
# Save the result of the feature selection
importance_df.to_csv(GB_working_dir + "feature_selection_results.csv", index = False)



