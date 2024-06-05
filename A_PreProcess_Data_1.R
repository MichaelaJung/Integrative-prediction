################################################################################
###                                                                          ###
### Preprocess data: Generating geno matrices, PCs, soil and weather         ###
### variables in a suitable format for the rest of the pipeline              ###                                                  
###                                                                          ###
### by Carles Quesada-Traver, 2024                                           ###
###                                                                          ###
################################################################################

# Usage: Run the script (i.e. Rscript A_PreProcess_Data_1.R)

# Clean the environment
rm(list = ls())

# Load required packages
library(SeqArray)                                                              
library(SNPRelate)                                                              
library(tidyverse)                                                              
library(GENESIS)
library(lme4)
library(readxl)


# Define variables to apply a filter of 97.5% LD and 2%MAF
LD  <- 0.975
MAF <- 0.02

# Define paths for inputs in plink format and output in gds format
bed <- 'Input/SNPs_final_2022.bed'
bim <- 'Input/SNPs_final_2022.bim'
fam <- 'Input/SNPs_final_2022.fam'
gds <- 'Output/Intermediate/SNPs_final_2022.gds'

# Transform to gds format for efficient data handling
seqBED2GDS(bed.fn = bed,
           fam.fn = fam,
           bim.fn = bim, 
           out.gdsfn = gds)

# Open a connection to the gds file
gds <- seqOpen(gds.fn = gds)

# Prune and filter the markers according to the above defined thresholds
snpset <- snpgdsLDpruning(gds, ld.threshold = LD,
                          autosome.only = FALSE, maf = MAF)
snpset <- unlist(unname(snpset))
seqSetFilter(gds, variant.id = snpset)

# Calculate the KING-robust estimate
ibd.robust <- snpgdsIBDKING(gds, snp.id = snpset, 
                            type = 'KING-robust', num.thread = 8)
KINGmat <- ibd.robust$kinship
colnames(KINGmat) <- rownames(KINGmat) <- ibd.robust$sample.id

# Calculate PCAir and save PC output for further processing
mypcair <- pcair(gds, kinobj = KINGmat, divobj = KINGmat, 
                 snp.include = snpset, eigen.cnt = 0)
write.csv(mypcair$vectors, quote = F, 
          file = 'Output/Intermediate/PCAir_PCs.csv')

# Create files to run GradientBoosting for feature selection using Python
# Prepare a common geno file (genotypes x snps)

# Extract genotype matrix from gds file
geno <- snpgdsGetGeno(gds, snp.id=snpset, snpfirstdim=FALSE)
geno <-  2 - geno
# Name the columns and rows according to the annotations
snps <- seqGetData(gds, var.name ="annotation/id")
genotypes <- seqGetData(gds, var.name ="sample.id")
colnames(geno) <- snps
rownames(geno) <- paste("G", genotypes, sep = "_")
# Save the geno file
write.csv(geno, file ="Output/Intermediate/GB_feature_selection/all.geno", quote =F, row.names = T) 

# Create a .csv for each combination of trait and split (n = 275)
# Load the phenotypic data
data <- as.data.frame(read_xlsx("Output/Pheno_processed.xlsx"))
data <- pivot_wider(data, names_from = Trait, values_from = Value)
# Save pivoted phenotypes for easier processing in python
write.csv(data, "Output/Intermediate/Pheno_processed_pivoted.csv", quote = F, row.names = F)
# Iterate over the 11 different traits
traits <- colnames(data)[3:13]
for (trait in traits) {
  # Filter for NA and order to match CV template
  data_trait <- na.omit(data[,c("Envir", "Genotype", trait)])
  data_trait <- data_trait[order(data_trait$Envir,data_trait$Genotype),]
  CV_strategy <- read.csv(paste0("Input/CV_Strategy/", trait,"_CV.csv"))
  CV_strategy <- CV_strategy[order(CV_strategy$Envir,CV_strategy$Genotype),]
  # Create a directory to save the 25 files if it does not exist
  save_dir = paste0("Output/Intermediate/GB_feature_selection/",trait)
  ifelse(!dir.exists(save_dir), dir.create(save_dir, showWarnings = F), F)
  # Iterate over splits
  splits <- colnames(CV_strategy)[3:27]
  for (split in splits) {
    # Filter for training split (0)
    pheno_train_split <-   data_trait[CV_strategy[split]== 0,]
    # Calculate BLUPs
    formula <- paste0(trait, " ~ Envir + (1|Genotype)")
    lmod <- lmer(formula, data = pheno_train_split)
    Y <- as.data.frame(ranef(lmod)$Genotype)
    Y <- cbind(rownames(Y), Y)
    colnames(Y) <- c("Genotype", trait)
    Y$Genotype <- paste("G", Y$Genotype, sep = "_")
    # Scale the data and save it for further processing
    Y[[trait]] <- scale(Y[[trait]])
    write.csv(Y,file =paste0(save_dir,"/", split, ".csv"),quote=F,row.names=F)
  } 
}

# Save .Rdata daily weather to .csv for python processing
load("Output/Weather_daily.Rdata")

# Remove some redundant information
d <- d[,1:5]
write.csv(d,"Output/Intermediate/Weather_daily.csv", quote = F, row.names = F)

# Save soil variables to .csv for python processing
load("Output/W.Rdata")
soil <-  W[,-c(9,10,23,24,25,26)]
soil <- cbind(Envir = rownames(soil), soil)
write.csv(soil, "Output/Intermediate/Soil.csv", quote = F, row.names = F)
