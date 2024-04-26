################################################################################
###                                                                          ###
### Genomic prediction integrating additive and dominance random effects     ###
### and inbreeding as fixed effect                                           ###
###                                                                          ###
### by Michaela Jung, 2024                                                   ###
###                                                                          ###
################################################################################

#############
### Setup ###
#############

# load libraries
library(snpStats)
library(readxl)
library(EnvRtype)
library(BGGE)

# make directory for the output
if (dir.exists("./Output/out_predictions")==FALSE) dir.create("./Output/out_predictions")


###############################
### Load data and functions ###
###############################

# source functions for the estimation of relationship matrices
source("Relationship_matrices_A_D.R")

# load genomic data
geno <- read.plink("./Input/SNPs_final_2022.bed",
                   "./Input/SNPs_final_2022.bim",
                   "./Input/SNPs_final_2022.fam")

# load phenotypic data
pheno <- as.data.frame(read_xlsx("./Output/Pheno_processed.xlsx"))

# # original settings (server)
# args <- commandArgs(trailingOnly = TRUE) # command line arguments
# eval(parse(text=args[1]))                # command line arguments
# eval(parse(text=args[2]))                # command line arguments
# eval(parse(text=args[3]))                # command line arguments
# nfolds <- 5       # number of folds in cross validation
# nIter <- 10000    # numer of iterations
# burnIn <- 1000    # burn-in
# nThin <- 3        # thinning
# testsetup <- "no" # run test setup following example settings (yes/no)

# example settings
var <- 1             # trait no. 1
run <- 1             # cross-validation repetition no. 1
mymodel <- 2         # model to run
nfolds <- 5          # number of folds in cross-validation 
nIter <- 50          # number of iterations
burnIn <- 10         # burn-in
nThin <- 2           # thinning
X_subs <- 10         # number of markers to subset
testsetup <- "yes"   # run test setup following example settings (yes/no)


#########################
### Format phenotypes ###
#########################

# subset for one trait
alltraits <- c("Harvest_date","Flowering_intensity","Fruit_weight",
               "Fruit_number","Fruit_weight_single","Flowering_begin",
               "Color_over","Russet_freq_all",
               "Acidity","Sugar","Firmness")
mytrait <- alltraits[var]
pheno <- pheno[which(pheno$Trait %in% mytrait),c("Envir","Genotype","Value")]
colnames(pheno)[3] <- "Trait"

# test setup - subset data for fewer genotypes
if (testsetup=="yes") { 
  pheno <- pheno[which(pheno$Genotype %in% unique(pheno$Genotype)[1:50]),]
  }

# format columns for BGGE
# the first column is a factor for environments
# the second column is a factor identifying genotypes
# the third column contains the trait
pheno$Envir <- as.factor(pheno$Envir)
pheno$Genotype <- as.factor(pheno$Genotype)


#######################
### Format SNP data ###
#######################

# obtain genotype values
X <- geno$genotypes@.Data

# test setup - subset matrix for fewer markers
if (testsetup=="yes") X <- X[,1:X_subs]

# numeric matrix of values 0, 1, 2
class(X) <- "numeric"
X <- X - 1

# subset for genotypes present in the phenotypic set
X <- X[which(row.names(X) %in% unique(as.character(pheno$Genotype))),]

# delete markers with no variation
v <- apply(X, 2, var)
if (length(which(v == 0)) > 0) X <- X[,-which(v == 0)]


##############################################
### Estimate genomic relationship matrices ###
##############################################

# measure time to estimate matrices
start_time1 <- Sys.time()

# additive genomic relationship matrix
A <- G.matrix(A.matrix(X))

# dominance genomic relationship matrix
if (mymodel %in% c(2, 4, 6)) D <- G.matrix(D.matrix(X))

# end the time measurement
end_time1 <- Sys.time()


#######################################
### Estimate inbreeding coefficient ###
#######################################

# set all homozygotes to zero
X2 <- X
X2[X2==2] <- 0

# inbreeding coefficient as the proportion of homozygotes
ic <- 1 - (rowSums(X2) / ncol(X2))

# merge inbreeding coefficient with the phenotypic data
# order genotypes as in the phenotypic data
ic_df <- data.frame(Genotype=names(ic), IC=ic)
pheno2 <- pheno
pheno2$id <- 1:nrow(pheno2)
pheno2 <- merge(pheno2, ic_df, all.x = TRUE)
pheno2 <- pheno2[order(pheno2$id), ]


##################
### Prediction ###
##################

# make a copy of the full phenotypic dataset
pheno_full <- pheno

# set random seed
set.seed(run)

# make a vector of all genotype names
genoID <- as.character(unique(pheno_full$Genotype))

# create folds
folds <- sample(cut(seq(1,length(genoID)),breaks=nfolds,labels=FALSE))

# initiate outputs
out1 <- NULL
out2 <- NULL
out3 <- NULL

### CROSS VALIDATION ###

# run a K-fold cross validation
# prediction of unknown genotypes across known environments
# all environments known
for (fold in 1:nfolds) {
  
  
  ### SETUP ###
  
  # test setup
  #fold <- 4
  
  # get the full set of phenotypic data
  pheno <- pheno_full
  
  # define the validation set - set part of genotypes to missing
  genoID_val <- genoID[which(folds == fold)]
  pheno[which(pheno$Genotype %in% genoID_val), "Trait"] <- NA
  
  # number of genotypes per environment
  ne <-  as.vector(table(pheno$Envir))
  
  
  ### FIXED EFFECTS ###
  
  # design matrix for the fixed effects of environments (all models)
  Z_E <- model.matrix(~0+Envir, pheno)
  
  # design matrix for the fixed effects of environments and inbreeding coefficient
  if(mymodel %in% c(5, 6)) Z_E <- model.matrix(~0+Envir+IC, pheno2)
  
  
  ### MODEL STRUCTURES ###
  
  # measure time to define model structures and fit the model
  start_time2 <- Sys.time()
  
  # no interaction effects
  if (mymodel %in% c(1, 5) ) {
    K_G <- list(A = A) # additive effects
    K <- get_kernel(K_G = K_G, K_E = NULL, env = 'Envir', gid='Genotype', 
                    y='Trait', data = pheno, model = 'MM')
  }
  
  if (mymodel %in% c(2, 6) ) {
    K_G <- list(A = A, D = D) # additive and dominance effects
    K <- get_kernel(K_G = K_G, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MM')
  }
  
  # with interaction effects
  if (mymodel == 3) {
    K_G <- list(A = A) # additive effects
    K <- get_kernel(K_G = K_G, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MDs')
  }
  
  if (mymodel == 4) {
    K_G <- list(A = A, D = D) # additive and dominance effects
    K <- get_kernel(K_G = K_G, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MDs')
  }

  
  ### PREDICTION ###
  
  # fit the model
  fit <- BGGE(y = pheno$Trait, K = K, XF= Z_E, ne = ne, ite = nIter,
              burn = burnIn, thin = nThin, verbose = FALSE)
  
  # end the time measurement
  end_time2 <- Sys.time()
  
  
  ### OUTPUT ###
  
  # sum both time measurements in seconds
  mytime1 <- as.numeric(difftime(end_time1, start_time1, units = "secs"))
  mytime2 <- as.numeric(difftime(end_time2, start_time2, units = "secs"))
  mytime <- mytime1 + mytime2
  
  # predicted values
  pheno$yHat <- as.vector(fit$yHat)
  
  # full set of phenotypic values
  pheno$yFull <- pheno_full$Trait
  
  # subset predicted and phenotypic values to the validation set
  pheno_val <- pheno[which(pheno$Genotype %in% genoID_val), ]
  
  # obtain summary statistics - variances and time
  if (mymodel %in% c(1, 5) ) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG_A", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel %in% c(2, 6) ) {
    out1 <- rbind(out1, c(fit$K$KG_G_A$varu, fit$K$KG_G_D$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG_A", "varG_D", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel == 3) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$K$KG_GE$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG_A", "varGE_A", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel == 4) {
    out1 <- rbind(out1, c(fit$K$KG_G_A$varu, fit$K$KG_G_D$varu,
                        fit$K$KG_GE_A$varu, fit$K$KG_GE_D$varu, fit$varE,
                        fold, mytime))
    colnames(out1) <- c("varG_A", "varG_D", "varGE_A", "varGE_D", "varRes",
                        "Fold", "Time")
  }
  
  # add identifiers for validation set and fold to the phenotypic data
  pheno$Validation <- NA
  pheno[which(pheno$Genotype %in% genoID_val), "Validation"] <- 1
  colnames(pheno)[3] <- "yTrain"
  pheno$Fold <- fold
  out2 <- rbind(out2, pheno)
  
  # estimate predictive ability as correlation separately for each environment
  mycor <- NULL
  for (i in unique(pheno_val$Envir)) {
    p <- pheno_val[which(pheno_val$Envir %in% i),]
    myc <- cor(p$yHat, p$yFull, use="pairwise.complete.obs")
    mycor <- rbind(mycor, c(myc, i))
  }
  mycor <- as.data.frame(mycor)
  mycor$Fold <- fold
  out3 <- rbind(out3, mycor)
  
}


##############
### Output ###
##############

# add identifiers for trait, model and repetition of cross-validation
out1 <- as.data.frame(out1)
out1$Trait <- mytrait
out2$Trait <- mytrait
out3$Trait <- mytrait

out1$Model <- mymodel
out2$Model <- mymodel
out3$Model <- mymodel

out1$Run <- run
out2$Run <- run
out3$Run <- run

# save output
write.table(out1, file = paste0("./Output/out_predictions/statistics_NOIA_", mytrait, 
                                "_model_", mymodel,
                                "_run", run, ".txt"),
            quote = FALSE, row.names = FALSE)

write.table(out2, file = paste0("./Output/out_predictions/predictions_NOIA_", mytrait,
                                "_model_", mymodel,
                                "_run_", run, ".txt"),
            quote = FALSE, row.names = FALSE)

colnames(out3)[1:2] <- c("r", "Envir")
write.table(out3, file = paste0("./Output/out_predictions/correlations_NOIA_", mytrait,
                                "_model_", mymodel,
                                "_run", run, ".txt"),
            quote = FALSE, row.names = FALSE)

