################################################################################
###                                                                          ###
### Relationship matrices for additive and dominance effects                 ###
###                                                                          ###
### Implementation of Álvarez-Castro & Carlborg (2007) Genetics              ###
### as detailed by Vitezica et al. (2017) Genetics                           ###
### https://doi.org/10.1534/genetics.116.199406                              ###
###                                                                          ###
### Adapted from Roth, Mary-Huard (2020)                                     ###
### https://doi.org/10.15454/ZGP766/0GACK9                                   ###
###                                                                          ###
### Adaptation of the G.matrix function from the R package snpReady          ###
### https://CRAN.R-project.org/package=snpReady                              ###
###                                                                          ###
### Adapted by Michaela Jung, 2024                                           ###
###                                                                          ###
################################################################################


### Matrix of coefficients for the additive contrasts

A.matrix <- function(M){ # M is a numeric matrix with genotypes (0, 1, 2)
  n <- nrow(M) # number of individuals
  Mdat <- lapply(as.data.frame(M), factor, levels=c("0","1","2"))
  p <- sapply(Mdat, function(x) table(x)) / n
  paa <- p[1,] # genotype 0
  pAa <- p[2,] # genotype 1
  pAA <- p[3,] # genotype 2
  S <- ((M==2)*1) * (-rep(-pAa-2*paa, each=n)) +
    ((M==1)*1) * (-rep(1-pAa-2*paa, each=n)) +
    ((M==0)*1) * (-rep(2-pAa-2*paa, each=n))
  return(S) 
}


### Matrix of coefficients for the dominant contrasts

D.matrix <- function(M){ # M is a numeric matrix with genotypes (0, 1, 2)
  n <- nrow(M) # number of individuals
  Mdat <- lapply(as.data.frame(M), factor, levels=c("0","1","2"))
  p <- sapply(Mdat, function(x) table(x)) / n
  paa <- p[1,] # genotype 0
  pAa <- p[2,] # genotype 1
  pAA <- p[3,] # genotype 2
  S <- ((M==2)*1) * (-rep(( (2*pAa*paa) / (pAA+paa-(pAA-paa)^2) ), each=n)) +
    ((M==1)*1) * rep(( (4*pAA*paa) / (pAA+paa-(pAA-paa)^2) ), each=n) +
    ((M==0)*1) * (-rep(( (2*pAA*pAa) / (pAA+paa-(pAA-paa)^2) ), each=n))
  S <- S[,which(!is.na(colSums(S)))] # removal of markers with NaN values
  return(S) 
}


### Genomic relationship matrix for the additive and dominance effects

G.matrix <- function(H) {
  n <- nrow(H) # number of individuals
  HHprime <- tcrossprod(H) 
  G <- HHprime*n / sum(diag(HHprime))
  return(G)
}

