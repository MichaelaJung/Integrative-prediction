################################################################################
###                                                                          ###
### Outlier detection methods for generalized lattices:                      ###
### A case study on the transition from ANOVA to REML                        ###
### Bernal-Vasquez et al. (2016) Theor Appl Genet                            ###
###                                                                          ###
### Adapted by Michaela Jung, 2024                                           ###
###                                                                          ###
################################################################################

outliers <- function(mydata, mytrait, mylm) {
  
  ##########################################################################
  ### METHOD 4: Bonferroni-Holm using re-scaled MAD for ####################
  ################# standardizing residuals (BH-MADR) ######################
  ##########################################################################
  
  ## Basic model
  
  # Analysis for original data
  #lmer.data=lmer(yield ~ gen + rep + (1|rep:block), data=dataset_name)
  lmer.data=mylm
  
  # re-scaled MAD
  resi=cbind(residuals(lmer.data, type = "response" ))
  median=median(resi)
  MAD=median((abs(resi - median)))
  re_MAD=MAD*1.4826
  # end
  
  # MAD standardized residuals
  res_MAD = resi /re_MAD
  # end
  
  # # Install package "multtest" from the bioconductor
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("multtest")
  library(multtest)
  
  # Calculate adjusted p-values
  rawp = 2 * (1 - pnorm(abs(res_MAD)))
  
  # Reduce the dataset in case of missing values.
  # These are not carried over to the rawp matrix
  #dataset_name.1= subset(dataset_name, dataset_name$yield!="NA")
  dataset_name.1= subset(mydata, mydata[,mytrait]!="NA")
  # Warning: Change variable that is under analysis, in this case yield
  
  # Combine the dataset, the residuals and the adjusted p-values
  rawp2 <- cbind(dataset_name.1, resi, res_MAD, rawp)
  
  # Produce a Bonferroni-Holm tests for the adjusted p-values
  # The output is a list
  res2<-mt.rawp2adjp(rawp,proc=c("Holm"))
  
  # Create vectors/matrices out of the list of the BH tests
  adjp = cbind(res2[[1]][,1])
  bholm = cbind(res2[[1]][,2])
  index = cbind(res2[[2]])
  
  # Condition to flag outliers according to the BH test
  out_flag = ifelse(bholm<0.05, "OUTLIER ", ".")
  
  # Create a matrix with all the output of the BH test
  bholm_test = cbind(adjp,bholm,index,out_flag)
  
  # Order the file by index
  bholm_test2 = bholm_test[order(index),]
  
  # Label colums
  names = c("rawp","bholm","index","out_flag")
  colnames(bholm_test2) <- names
  
  # Create a final file, with the data and the test and the labels for the outliers
  total.m4_data <- cbind(rawp2,bholm_test2)
  
  # Take a look at the outliers
  return(outliers_BHMAD <- total.m4_data[which(total.m4_data$out_flag!="."),])
  
  ##########################################################################

  
}

