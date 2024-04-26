################################################################################
###                                                                          ###
### Relationship matrices for genomic and enviromic effects                  ###
###                                                                          ###
### From Costa-Neto et al. (2021) Heredity                                   ###
### https://doi.org/10.1038/s41437-020-00353-1                               ###
### Code available at: https://github.com/gcostaneto/KernelMethods           ###
###                                                                          ###
### Adapted by Michaela Jung, 2024                                           ###
###                                                                          ###
################################################################################


### Relationship matrix based on G-BLUP approach

GB_Kernel <- function(X, is.center=FALSE) {
  if(isFALSE(is.center)) X = scale(x = X, center = T, scale = F)
  XXl <- X %*% t(X)
  K_G <- XXl/(sum(diag(XXl))/nrow(X)) + diag(1e-6, nrow(XXl))
  return(K_G)
}


### Relationship matrix based on Gaussian kernel

GK_Kernel <- function (X, h = NULL, prob=.5, y=NULL) {
  
  margh.fun <- function(theta, y, D, q, nu=0.0001, Sc=0.0001, nuh=NULL, Sch=NULL, prior=NULL) {
    
    h <- theta[1]
    phi <- theta[2]
    Kh <- exp(-h*D/q)
    eigenKh <- eigen(Kh)
    nr <- length(which(eigenKh$val> 1e-10))
    Uh <- eigenKh$vec[,1:nr]
    Sh <- eigenKh$val[1:nr]
    d <- t(Uh) %*% scale(y, scale=F)
    
    Iden <- -1/2*sum(log(1+phi*Sh)) - (nu+nr-1)/2*log(Sc+sum(d^2/(1+phi*Sh)))
    if(!is.null(prior)) lprior <- dgamma(h,nuh,Sch,log=T) else lprior <- 0
    
    Iden <- -(Iden+lprior)
    
    return(Iden)
  }
  
  GK <-list()
  for(i in 1:length(X)){
    d <- as.matrix(dist(X[[i]], upper = T, diag = T))^2
    q <- quantile(x=d,prob=prob)
    if(q == 0) q <- quantile(x=d,prob=.05)
    if (is.null(h)) h <- 1
    
    if(!is.null(y)){
      sol<-optim(c(1,3),margh.fun,y=y,D=d,q=q,method="L-BFGS-B",
                 lower=c(0.05,0.05),upper=c(6,30),prior="gamma",nuh=3,Sch=1.5)
      h<-sol$par[1]
    }
    
    GK[[i]] <- exp(-h * d/q)
  }
  names(GK) <- names(X)
  return(GK)
}


### Relationship matrix based on Deep kernel

get_GC1 <- function(M){ # Author: Germano Costa Neto
  AK1 <- list()
  for(i in 1:length(M)) AK1[[i]] <- GC1.fun(X = M[[i]])
  length(AK1)
  names(AK1) = names(M)
  return(AK1)
}

GC1.fun <- function(X){ # Author: Jaime Cuevas
  n <- nrow(X)
  cosalfa <- cor(t(X))
  angulo <- acos(cosalfa)
  mag <- sqrt(apply(X,1,function(x) crossprod(x)))
  sxy <- tcrossprod(mag)
  GC1 <- (1/pi)*sxy*(sin(angulo)+(pi*matrix(1,n,n)-angulo)*cosalfa)
  GC1 <- GC1/median(GC1)
  colnames(GC1) <- rownames(X)
  rownames(GC1) <- rownames(X)
  return(GC1)
}  

