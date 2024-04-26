################################################################################
###                                                                          ###
### Phenotypic data preprocessing                                            ###
###                                                                          ###
### by Michaela Jung, 2024                                                   ###
###                                                                          ###
################################################################################

#############
### Setup ###
#############

# load libraries
library(readxl)
library(lme4)
library(snpStats)
library(lattice)
library(openxlsx)
library(SpATS)


#################
### Load data ###
#################

# load phenotypic data
pheno_raw <- as.data.frame(read_xlsx("./Input/Pheno_raw.xlsx"))

# load genomic data
geno <- read.plink("./Input/SNPs_final_2022.bed",
                   "./Input/SNPs_final_2022.bim",
                   "./Input/SNPs_final_2022.fam")

# traits to analyze
traits <- c("Harvest_date","Fruit_weight","Fruit_number","Fruit_weight_single",
            "Color_over","Russet_freq_all",
            "Flowering_intensity","Flowering_begin",
            "Acidity","Sugar","Firmness")


###########################
### Data transformation ###
###########################

# log-transformation of yield-related traits
pheno_raw$Fruit_weight <- log(pheno_raw$Fruit_weight)
pheno_raw$Fruit_number <- log(pheno_raw$Fruit_number)


####################
### Heritability ###
####################

# initialize output
her <- NULL

# loop over traits
for (i in traits) {

  # subset for one trait
  d <- pheno_raw[,c("Envir","Genotype",i)]
  colnames(d)[3] <- "Variable"

  # filter for genotypes with genomic data
  d <- d[which(d$Genotype %in% geno$fam$member),]

  # loop over environments
  for (j in unique(d$Envir)) {

    # subset for one environment
    d_subs <- subset(d, Envir %in% j)

    # remove environments with too few or no observations
    # skip to next iteration if TRUE
    if (sum(is.na(d_subs$Variable)) > 1100) {
      next
    }

    # remove missing values
    if (any(is.na(d_subs$Variable))) {
      d_subs <- d_subs[-which(is.na(d_subs$Variable)),]
    }

    # variables as factors
    d_subs$Genotype <- as.factor(as.character(d_subs$Genotype))
    d_subs$Envir <- as.factor(as.character(d_subs$Envir))

    # fit mixed-effects model
    mylm <- lmer(Variable ~ (1|Genotype), data = d_subs)

    # extract variance components
    vc <- as.data.frame(VarCorr(mylm))

    # genotypic variance
    genov <- vc$vcov[vc$grp=="Genotype"]

    # error variance
    errorv <- sigma(mylm)^2

    # mean number of replicates
    mean_reps <- mean(replications(d_subs)$Genotype)

    # masking variance
    maskv <- errorv/mean_reps

    # heritability
    H2 <- genov / (genov + maskv)

    # save output
    her <- rbind(her,data.frame(i, H2=H2, Envir=j))

  }

}

# save H2
dir.create(file.path("./Output"))
save(her, file = "./Output/heritability.RData")


#############################
### Spatial heterogeneity ###
#############################

# initialize output
trees.all <- NULL
means.all <- NULL

# variables to keep in dataframe
ids <- c("Row","Position","Genotype","Management","Checks","Year","Country","Envir")

# loop over runs
for (envir in unique(pheno_raw$Envir)) {
  
  # subset data
  data <- pheno_raw[which(pheno_raw$Envir %in% envir),]
  
  # input variables for the models
  data$column <- as.numeric(data$Row)
  data$row <- as.numeric(data$Position)
  data$R <- as.factor(data$row)
  data$C <- as.factor(data$column)
  data$Genotype <- as.factor(data$Genotype)
  
  # initialize output
  trees <- data[,ids]
  trees <- trees[order(trees$Row, trees$Position),]
  means <- data.frame(Genotype=as.character(unique(data$Genotype)),
                      stringsAsFactors=FALSE)
  
  # loop over traits
  for (mytrait in traits) {
    
    # remove environments with too few or no observations
    # skip to next iteration if TRUE
    if (sum(is.na(data[,mytrait])) > 1100) {
      next
    }
    
    # subset data frame for one trait
    d <- data[,c("column", "row", "C", "R", ids, mytrait)]
    colnames(d)[ncol(d)] <- "Trait"
    
    # run the spatial analysis with splines
    spatial.SpATS <- as.formula(~ PSANOVA(column, row))
    m <- SpATS(response = "Trait",
               spatial = spatial.SpATS,
               genotype = "Genotype", genotype.as.random = TRUE, fixed = NULL,
               random = ~ C + R, data = d,
               control =  controlSpATS())
    
    # prediction from a SpATS object
    pred.m <- predict(m, which = "Genotype")
    pred.m <- pred.m[,c("Genotype","predicted.values")]
    
    # adjusted phenotypic values of each genotype
    means <- merge(means, pred.m, all.x=TRUE)
    colnames(means)[ncol(means)] <- mytrait
    
    # residuals
    d1 <- m$data
    d1$weights <- NULL
    colnames(d1)[ncol(d1)] <- mytrait
    d1$residuals <- m$residuals
    
    # merge predicted values with residuals
    d1$Genotype <- as.character(d1$Genotype)
    pred.m$Genotype <- as.character(pred.m$Genotype)
    d2 <- merge(d1, pred.m, all.x = TRUE)
    
    # delete predicted values that were originally missing
    d2[which(is.na(d2[,mytrait])),"predicted.values"] <- NA
    
    # adjusted phenotypic values of each tree
    d2$adjusted <- d2$predicted.values + d2$residuals
    
    # order data by row and position of trees
    d2 <- d2[order(d2$Row, d2$Position),]
    
    # save results into a common data frame
    trees <- cbind(trees, d2$adjusted)
    colnames(trees)[ncol(trees)] <- mytrait
    
  }
  
  # save output per tree
  trees.all <- dplyr::bind_rows(trees.all, trees)
  
  # save output per genotype
  means <- cbind(rep(envir, nrow(means)), means)
  colnames(means)[1] <- "Envir"
  means.all <- dplyr::bind_rows(means.all, means)
  
}

# save
save(trees.all, file = "./Output/adjusted_values_trees.RData")
save(means.all, file = "./Output/adjusted_values_genotype.RData")


######################
### Data filtering ###
######################

# load previously produced output for spatial heterogeneity
# load("./Output/adjusted_values_genotype.RData")

# melt dataframe
pheno <- means.all
pheno <- reshape2::melt(pheno, id.vars=c("Envir","Genotype"))
colnames(pheno)[3:4] <- c("Trait", "Value")
pheno$Trait <- as.character(pheno$Trait)

# remove observations for genotypes without genomic data
pheno <- pheno[which(pheno$Genotype %in% geno$fam$pedigree),]

# remove missing values
pheno <- pheno[-which(is.na(pheno$Value)),]


### remove trait-environment combinations with low heritability (H2<0.1) ###

# low heritability
her <- her[which(her$H2 < 0.1),]
her

# for the environment with low heritability for floral emergence,
# estimate the date when 90% of genotypes started flowering
round(quantile(
  pheno[which(pheno$Envir %in% "ESP.2020" & pheno$Trait %in% "Flowering_begin"),
        "Value"],
  probs = 0.9, na.rm=T))

# filter for heritability
for (i in 1:nrow(her)) {
  trait <- her[i,1]
  envir <- her[i,3]
  pheno <- pheno[-which(pheno$Envir %in% envir & pheno$Trait %in% trait),]
}


################
### Outliers ###
################

# source code for outlier detection
source("Outliers.R")

# make directory for graphical output
dir.create("./Output/out_outliers")

# start output
pheno_out <- NULL

# loop over traits
for (mytrait in traits) {
  
  # subset for a trait
  p <- pheno[which(pheno$Trait %in% mytrait),]
  
  # edit rownames
  rownames(p) <- 1:nrow(p)
  
  # format predictors as factors
  p$Genotype <- as.factor(p$Genotype)
  p$Envir <- as.factor(p$Envir)
  
  # fit the mixed-effects model
  mylm <- lmer(Value ~  Envir + (1|Genotype), data = p)
  
  # call function for outlier identification
  outl <- outliers(p, "Value", mylm)
  
  # formatting of output
  d <- p
  d$Predicted <- NA
  d[names(predict(mylm)),"Predicted"] <- predict(mylm)
  d$Outlier <- NA
  d[rownames(outl),"Outlier"] <- d[rownames(outl),"Value"]
  
  # plot
  jpeg(filename=paste0("./Output/out_outliers/Outliers_",mytrait,".jpg"),
       width = 15, height = 15, units = "cm", res=300)
  print(plot(d$Value ~ d$Predicted))
  points(d$Outlier ~ d$Predicted, col="red")
  dev.off()
  
  # remove outliers
  if (nrow(outl) > 0) p <- p[-as.numeric(rownames(outl)),]
  
  # save removed outliers
  write.table(outl, paste0("./Output/out_outliers/Outliers_",mytrait,".txt"),
              col.names = TRUE, row.names=FALSE, quote = FALSE)
  
  # save subset of phenotypes
  pheno_out <- rbind(pheno_out, p)
  
}

# save output
write.xlsx(pheno_out, "./Output/Pheno_processed.xlsx")


########################
### Phenotypic model ###
########################

# load previously produced output for spatial heterogeneity
# load("./Output/adjusted_values_trees.RData")

# melt dataframe
trees.all <- trees.all[,c("Envir","Genotype",traits)]
trees.all <- reshape2::melt(trees.all, id.vars=c("Envir","Genotype"))
colnames(trees.all)[3:4] <- c("Trait", "Value")
trees.all$Trait <- as.character(trees.all$Trait)

# remove observations for genotypes without genomic data
trees.all <- trees.all[which(trees.all$Genotype %in% geno$fam$pedigree),]

# remove missing values
trees.all <- trees.all[-which(is.na(trees.all$Value)),]

# remove trait-environment combinations with low heritability (H2<0.1)
for (i in 1:nrow(her)) {
  trait <- her[i,1]
  envir <- her[i,3]
  trees.all <- trees.all[-which(trees.all$Envir %in% envir & trees.all$Trait %in% trait),]
}

# start output
var <- NULL

# loop over traits
for (mytrait in traits) {

  # subset for a trait
  p <- trees.all[which(trees.all$Trait %in% mytrait),]
  
  # format predictors as factors
  p$Genotype <- as.factor(p$Genotype)
  p$Envir <- as.factor(p$Envir)

  # fit the mixed-effects model
  mylm <- lmer(Value ~  Envir + (1|Genotype) + (1|Genotype:Envir), data = p)

  # extract variance of the random effects
  vars_rand <- data.frame(VarCorr(mylm))

  # proportion of variance of the random effects
  vars_rand$Prop <- NA
  for (i in 1:nrow(vars_rand)) vars_rand[i,"Prop"] <- (vars_rand$vcov[i] ) / sum(vars_rand$vcov)

  # save output into a common data frame
  vars_rand$Trait <- mytrait
  var <- rbind(var, vars_rand)
  
}

# save
save(var, file="./Output/Phenotypic_model_variance.RData")

