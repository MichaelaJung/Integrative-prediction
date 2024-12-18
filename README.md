# Integrative-prediction

## General description

This respository contains the supporting code for the article "Integrative multi-environmental genomic prediction in apple" by Jung, Quesada-Traver, Roth, et al. (2024) _Horticulture Research_. It consists of two sets of scripts to be run sequentially (1-3; A-G) as well as scripts containing helper functions. Both sets perform genomic predictions of eleven different traits in apple using two different approaches: statistical methods (scripts 1-3) and deep learning (scripts A-G).

## Input data

To ensure code functionality, the phenotypic, enviromic, and imputed genomic data formatted as input files can be downloaded into the folder Input from Zenodo at https://doi.org/10.5281/zenodo.14191209. The folder Input additionally contains the Supplementary table 3 from Jung et al. (2022) _Horticulture Research_ listing SNPs associated with the studied traits. The subfolder CV1_Strategy contains the applied cross-validation strategy assigning individual genotypes into folds (i.e., splits).

## Statistical methods for genomic prediction

The outcomes of the tested statistical genomic prediction models can be reproduced by sequentially running the scripts indexed from 1 to 3. The first two scripts perform operations to process data into the format required for model training. The third script performs genomic predictions integrating random effects for the standard genomic, additive genomic, dominance genomic, and enviromic relationship matrix and the fixed effect of inbreeding.

## Deep learning method for genomic prediction

The deep learning implementation can be reproduced by sequentially running the scripts alphabetically indexed from A to G. Briefly, the first scripts (A-D) perform feature selection and data processing operations to accomodate it into the appropriate format for model training. Script E performs the crossvalidation for a given trait. Finally, scripts E-F compute and save in a comprehensive format interpretation metrics based on Shapley additive explanations (SHAP). Information on how to run each script as well as a short description of inputs and outputs is provided in their headers.

A representation of the deep model architecture generated with Netron (Netron: Visualizer for neural network, deep learning and machine learning models. https://www.lutzroeder.com/ai):

![Deep_learning_architechture](images/model_architecture_horiz.svg)
