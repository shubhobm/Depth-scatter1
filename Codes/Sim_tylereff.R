## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(parallel)
library(doSNOW)

## setup 1: Sigma = diag(1,1)
set.seed(02032015)

n = 50
Sig = diag(c(2,1))
simiter = 100
err = matrix(0, nrow=simiter, ncol=2)

for(i in 1:simiter){
  iX = matrix(rnorm(2*n), ncol=2) %*% sqrt(Sig)
  
  err[i,1] = norm(Sig - TylerSigma(iX), type="F")
  err[i,2] = norm(Sig - TylerSigma(iX, depth=T), type="F")
}

colSums(err)/simiter
