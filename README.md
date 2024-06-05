# Integrative-prediction

## General description

This respository contains the supporting code for the article "Integrative multi-environmental genomic prediction in apple". It consists of two sets of scripts to be run sequentially (1-4; A-G) as well as scripts containing helper functions. Both sets perform genomic predictions of eleven different traits in apple using two different approaches: classical statistical methods and deep learning respectively. The required files to run the pipelines from the begining are located in the 'Input' deirectory. However, the output and intermediate output is excluded and can be downloaded here (add hyperlink to zenodo repository with the Output files to download).

## Statistical methods for genomic prediction

## Deep learning method for genomic prediction

The deep learning implementation can be reproduced by sequentially running the scripts alphabetically indexed from A to G. Briefly, the first scripts (A-D) perform feature selection and data processing operations to accomodate it into the appropiate format for model training. Script E performs the crossvalidation for a given trait. Finally, scripts E-F compute and save in a comprehensive format interpretation metrics based on Shapely additive explanations (SHAP). Information on how to run each script as well as a short description of inputs and outputs is provided in their headers. 

A representation of the deep model architecture generated with Netron (Netron: Visualizer for neural network, deep learning and machine learning models. https://www.lutzroeder.com/ai):

![Deep_learning_architechture](images/model_architecture_horiz.svg)