## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')

library(fda.usc)
library(MethylCapSig)
library(parallel)

get.power = function(n,p,nrep,mu,cov,rho=0.5){
  
  ## construct output matrix
  out.mat = matrix(0, nrow=length(n)*length(p), ncol=5)
  out.mat[,1] = rep(p,rep(length(n),length(p)))
  out.mat[,2] = rep(n, length(p))
  
  ## get rows of output matrix, each by nrep simulations
  loopfun = function(j){
    
    ## Get null distributions
    is.reject = rep(0,3)
    
    # generate data
    X1 = matrix(rnorm(ni*pi), ncol=pi) %*% t(Sigma.sqrt)
    X2 = 3 * matrix(rnorm(ni*pi), ncol=pi) %*% t(Sigma.sqrt)
    ind = rbinom(ni, 1, 0.5)
    X = mu + ind*X1 + (1-ind)*X2
    
    # Chen-Qin test
    is.reject[1] = (cqtest(X)[2] < 0.05)
    
    # Sign-based test
    uX = X / sqrt(rowSums(X^2))
    is.reject[2] = (cqtest(uX)[2] < 0.05)
    
    # Depth-rank based test
    dX = mdepth.TD(X,X)$dep
    Xrank = uX * (max(dX)-dX)
    #Xrank = uX * dX
    is.reject[3] = (cqtest(Xrank)[2] < 0.05)
    
    is.reject
  }
  
  for(irep in 1:nrow(out.mat)){
    ni = out.mat[irep,2]
    pi = out.mat[irep,1]
    Sigma = diag(rep(1,pi))
    if(cov=="AR1"){
      for(i in 2:pi){
        for(j in 1:i-1){
          Sigma[i,j] = rho^abs(i-j)
          Sigma[j,i] = Sigma[i,j]
        }
      }
    }
    else if(cov=="compound"){
      for(i in 2:pi){
        for(j in 1:i-1){
          Sigma[i,j] = rho
          Sigma[j,i] = Sigma[i,j]
        }
      }
    }
    else{
      D = diag(2 + (pi-(1:pi)+1)/pi)
      R = diag(pi)
      for(i in 2:pi){
        for(j in 1:i-1){
            R[i,j] = -1^(i+j)*0.2^(abs(i-j)^0.1)
            R[j,i] = R[i,j]
        }
      }
      Sigma = D %*% R %*% D
    }
    eo = eigen(Sigma)
    Sigma.sqrt = eo$vectors %*% diag(sqrt(eo$values))
    
    power.vec = lapply(1:nrep, loopfun)
#    power.vec = mclapply(1:nrep, loopfun, mc.cores=16)
    power.vec = apply(matrix(unlist(power.vec), ncol=3, byrow=T), 2, mean)
    out.mat[irep,3:5] = power.vec
  }
  
  out.frame = data.frame(out.mat)
  names(out.frame) = c("p","n", "T2","Sign","Rank")
  out.frame
}

n = c(20,50)
p = c(500,1e3)
nrep = 1e2

get.power(n, p, nrep, mu=0.15, cov="AR1", rho=0.8)
#get.power(n, p, nrep, mu=0, cov="compound", rho=0.2)
#get.power(n, p, nrep, mu=0.2, cov="srivatsa")

