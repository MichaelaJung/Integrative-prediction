################################################################################
###                                                                          ###
### Genomic prediction integrating genomic and enviromic random effects      ###
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
source("Relationship_matrices_G_W.R")

# load genomic data
geno <- read.plink("./Input/SNPs_final_2022.bed",
                   "./Input/SNPs_final_2022.bim",
                   "./Input/SNPs_final_2022.fam")

# load phenotypic data
pheno <- as.data.frame(read_xlsx("./Output/Pheno_processed.xlsx"))

# load environmental covariables
load("./Output/W.Rdata")

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


#################################
### Format environmental data ###
#################################

# subset for environments present in the phenotypic dataset
envs <- unique(pheno$Envir)
W <- W[which(rownames(W) %in% envs),]

# scale and center the data
W <- scale(W, center=TRUE, scale=TRUE)

# remove columns with all values missing 
W <- W[, colSums(is.na(W)) != nrow(W)]


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

# scale and center
X <- scale(X, center = TRUE, scale = TRUE)


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
folds <- sample(cut(seq(1,length(genoID)), breaks=nfolds, labels=FALSE))

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
  
  # design matrix for the fixed effects of environments
  Z_E <- model.matrix(~0+Envir,pheno)
  
  
  ### MODEL STRUCTURES ###
  
  # measure time to define model structures and fit the model
  start_time <- Sys.time()
  
  # relationship matrices based on G-BLUP
  if (mymodel == 1) {
    K_A <- GB_Kernel(X = X, is.center = TRUE)
    K_A <- list(A = K_A)
    K <- get_kernel(K_G = K_A, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MM')
  }
  if (mymodel == 2) {
    K_A <- GB_Kernel(X = X, is.center = TRUE)
    K_A <- list(A = K_A)
    K <- get_kernel(K_G = K_A, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MDs')
  }
  if (mymodel == 3) {
    K_A <- GB_Kernel(X = X, is.center = TRUE)
    K_A <- list(A = K_A)
    K_W <- list(W=env_kernel(env.data = as.matrix(W))[[2]])
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'EMM')
  }
  if (mymodel == 4) {
    K_A <- GB_Kernel(X = X, is.center = TRUE)
    K_A <- list(A = K_A)
    K_W <- list(W=env_kernel(env.data = as.matrix(W))[[2]])
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'RNMM')
  }
  if (mymodel == 5) {
    K_A <- GB_Kernel(X = X, is.center = TRUE)
    K_A <- list(A = K_A)
    K_W <- list(W=env_kernel(env.data = as.matrix(W))[[2]])
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'RNMDs')
  }
  
  # relationship matrices based on Gaussian kernel
  if (mymodel == 6) {
    K_G   <- GK_Kernel(X = list(A=X))
    K_A   <- list(A = K_G$A)
    K <- get_kernel(K_G = K_A, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MM')
  }
  if (mymodel == 7) {
    K_G   <- GK_Kernel(X = list(A=X))
    K_A   <- list(A = K_G$A)
    K <- get_kernel(K_G = K_A, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MDs')
  }
  if (mymodel == 8) {
    K_G   <- GK_Kernel(X = list(A=X))
    K_A   <- list(A = K_G$A)
    K_W <- list(W=env_kernel(env.data = as.matrix(W), gaussian = TRUE)[[2]])
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'EMM')
  }
  if (mymodel == 9) {
    K_G   <- GK_Kernel(X = list(A=X))
    K_A   <- list(A = K_G$A)
    K_W <- list(W=env_kernel(env.data = as.matrix(W), gaussian = TRUE)[[2]])
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'RNMM')
  }
  if (mymodel == 10) {
    K_G   <- GK_Kernel(X = list(A=X))
    K_A   <- list(A = K_G$A)
    K_W <- list(W=env_kernel(env.data = as.matrix(W), gaussian = TRUE)[[2]])
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'RNMDs')
  }
  
  # relationship matrices based on Deep kernel
  if (mymodel == 11) {
    K_A  <- get_GC1(M = list(A=X))
    K <- get_kernel(K_G = K_A, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MM')
  }
  if (mymodel == 12) {
    K_A  <- get_GC1(M = list(A=X))
    K <- get_kernel(K_G = K_A, K_E = NULL, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'MDs')
  }
  if (mymodel == 13) {
    K_A  <- get_GC1(M = list(A=X))
    K_W  <- get_GC1(M = list(W=W))
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'EMM')
  }
  if (mymodel == 14) {
    K_A  <- get_GC1(M = list(A=X))
    K_W  <- get_GC1(M = list(W=W))
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'RNMM')
  }
  if (mymodel == 15) {
    K_A  <- get_GC1(M = list(A=X))
    K_W  <- get_GC1(M = list(W=W))
    K <- get_kernel(K_G = K_A, K_E = K_W, env = 'Envir', gid='Genotype',
                    y='Trait', data = pheno, model = 'RNMDs')
  }
  
  
  ### PREDICTION ###
  
  # fit the model
  fit <- BGGE(y = pheno$Trait, K = K, XF= Z_E, ne = ne, ite = nIter,
              burn = burnIn, thin = nThin, verbose = FALSE)
  
  # end time measurement
  end_time <- Sys.time()

  
  ### OUTPUT ###
  
  # estimate time in seconds
  mytime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # predicted values
  pheno$yHat <- as.vector(fit$yHat)
  
  # full set of phenotypic values
  pheno$yFull <- pheno_full$Trait
  
  # subset predicted and phenotypic values to the validation set
  pheno_val <- pheno[which(pheno$Genotype %in% genoID_val), ]
  
  # obtain summary statistics - variances and time
  if (mymodel %in% c(1,6,11)) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel %in% c(2,7,12)) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$K$KG_GE$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG", "var_GE", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel %in% c(3,8,13)) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$K$KE_W$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG", "var_W", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel %in% c(4,9,14)) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$K$KE_W$varu, fit$K$KGE_AW$varu,
                          fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG", "var_W","var_GW", "varRes",
                        "Fold", "Time")
  }
  
  if (mymodel %in% c(5,10,15)) {
    out1 <- rbind(out1, c(fit$K$KG_G$varu, fit$K$KG_GE$varu, fit$K$KE_W$varu,
                          fit$K$KGE_AW$varu, fit$varE,
                          fold, mytime))
    colnames(out1) <- c("varG", "var_GE","var_W","var_GW", "varRes",
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
write.table(out1, file = paste0("./Output/out_predictions/statistics_G_W_", mytrait, 
                                "_model_", mymodel,
                                "_run", run, ".txt"),
            quote = FALSE, row.names = FALSE)

write.table(out2, file = paste0("./Output/out_predictions/predictions_G_W_", mytrait,
                                "_model_", mymodel,
                                "_run_", run, ".txt"),
            quote = FALSE, row.names = FALSE)

colnames(out3)[1:2] <- c("r", "Envir")
write.table(out3, file = paste0("./Output/out_predictions/correlations_G_W_", mytrait,
                                "_model_", mymodel,
                                "_run", run, ".txt"),
            quote = FALSE, row.names = FALSE)

