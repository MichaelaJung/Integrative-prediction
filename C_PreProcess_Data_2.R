################################################################################
###                                                                          ###
### Preprocess data nr. 2: Processing the output from the GB and generating  ###
### suitable inputs for the rest of the pipeline                             ###
###                                                                          ###
### by Carles Quesada-Traver, 2024                                           ###
###                                                                          ###
################################################################################

# Usage: Run the script (i.e. "Rscript C_PreProcess_Data_2.R")

# Clean the environment
rm(list = ls())

# Load required packages
library(SeqArray)                                                              
library(SNPRelate)                                                              
library(tidyverse)
library(readxl)



importance_all <- read.csv("Output/Intermediate/GB_feature_selection/feature_selection_results.csv")

# Load GWAS SNPs that have been previously linked to studied traits
gwas_snps <- read_xls("Input/SupTable3_SNPS_GWAS.xls")[,1:2]

# Rename traits to match the trait names in the analysis
traits <- unique(gwas_snps$Trait)
traits_ord <- traits[order(traits)]
renamed_traits_ord <- c("Flowering_begin", "Flowering_intensity", "Firmness", "Harvest_date", "Fruit_number",
                        "Color_over", "Russet_freq_all", "Fruit_weight_single", "Sugar", "Acidity", 
                        "Fruit_weight")
conversion_table <-  as.data.frame(cbind(traits_ord,renamed_traits_ord))
gwas_snps <- merge(gwas_snps,conversion_table, by.x = "Trait", by.y = "traits_ord", all.x = TRUE)
gwas_snps <- distinct(gwas_snps)
gwas_snps <- gwas_snps[,-1]

# Open connection to the gds file
gds_path <- 'Output/Intermediate/SNPs_final_2022.gds'

#Extract SNP IDs and Genotype IDs for filtering and annotating geno files per split
gds <- seqOpen(gds.fn = gds_path)
snp_idx_gds<- seqGetData(gds, var.name ="annotation/id")
genotypes <- seqGetData(gds, var.name ="sample.id")
seqClose(gds)

# Now iterating over traits and splits
traits <- unique(importance_all$Trait)
splits <- unique(importance_all$Split)
for (trait in traits) {
  # Store SNPs to add from GWAS for the given trait
  snps_gwas_add <- gwas_snps[gwas_snps$renamed_traits_ord == trait,]$SNP
  
  # Define a directory to save the geno .csv files
  save_dir = paste0("Output/Intermediate/geno_files/",trait)
  ifelse(!dir.exists(save_dir), dir.create(save_dir, showWarnings = F), F)
  for (split in splits) {
    gds <- seqOpen(gds.fn = gds_path)
    # Filter general importance df for trait and split of interest
    snps_important_df  <- importance_all[importance_all$Trait == trait & importance_all$Split == split,]
    
    #Order them according to decreasing importance
    snps_important_df <- snps_important_df[order(snps_important_df$importance, decreasing = TRUE),]
    
    # Select 1000 most important SNPs
    snps_important <- snps_important_df$SNP[1:1000]
    
    # Add the relevant GWAS snps for given trait
    snps_important_all <- unique(c(snps_important,snps_gwas_add))
    
    # Obtain the indices of the selected snps from the gds file
    snp.id <- (1:length(snp_idx_gds))[snp_idx_gds %in% snps_important_all]
    
    # Filter the gds file and generate the genotypic matrix
    geno_split <- snpgdsGetGeno(gds, snp.id=snp.id, snpfirstdim=FALSE)
    geno_split <- (2 - geno_split)/2
    
    # Add SNP and genotype names and save the matrix
    colnames(geno_split) <- snp_idx_gds[snp.id]
    rownames(geno_split) <- paste("G", genotypes, sep = "_")
    write.csv(geno_split, paste0(save_dir,"/geno_",split,".csv"), quote = F, row.names = T)
    seqClose(gds)
  }
}