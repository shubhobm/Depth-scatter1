## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')

library(fda.usc)
library(WCQ)
library(MethylCapSig)

level.vec = rep(0,3)
mu = .25
nrep = 1e2
n = 20
p = 1e3
rho = 0.8
Sigma = diag(rep(1,p))
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = rho^abs(i-j)
  }
}
eo = eigen(Sigma)
Sigma.sqrt = eo$vectors %*% diag(sqrt(eo$values))

level.vec = rep(0,3)
## Crossprod hiigh-dimensional test
pb = txtProgressBar(0,nrep)
for(rep in 1:nrep){
  
  # generate data
  #set.seed(1)
  X1 = matrix(rnorm(n*p), ncol=p) %*% t(Sigma.sqrt)
  X2 = 3 * matrix(rnorm(n*p), ncol=p) %*% t(Sigma.sqrt)
  X = mu + .9*X1 + .1*X2
#   uX = X / sqrt(rowSums(X^2))
#   dX = mdepth.TD(X,X)$dep
#   Xrank = uX * exp(-dX)
#   cqtest(X)
#   cqtest(uX)
#   cqtest(Xrank)
  
  # Chen-Qin test
  level.vec[1] = level.vec[1] + (cqtest(X)[2] < 0.05)
  
  # Sign-based test
  uX = X / sqrt(rowSums(X^2))
  level.vec[2] = level.vec[2] + (cqtest(uX)[2] < 0.05)
  
  # Depth-rank based test
  dX = mdepth.TD(X,X)$dep
  Xrank = uX * (max(dX)-dX)
  #Xrank = uX * dX
  level.vec[3] = level.vec[3] + (cqtest(Xrank)[2] < 0.05)
  setTxtProgressBar(pb,rep)
}
close(pb)

level.vec = level.vec/nrep
level.vec 

n=1e3
p=10
mu=0
rho = 0.8
Sigma = diag(rep(1,p))
for(i in 1:p){
  for(j in 1:p){
    Sigma[i,j] = rho^abs(i-j)
  }
}

X1 = my.mvrnorm(n, mu=rep(mu,p), Sig=Sigma)
X2 = my.mvrnorm(n, mu=rep(mu,p), Sig=9*Sigma)
X = .9*X1 + .1*X2

# bootstrap sign test
uX0 = (X-mu) / sqrt(rowSums((X-mu)^2))
uX0bar = apply(uX0,2,mean)

nboot=1e3
Xb = X - matrix(apply(X,2,mean), ncol=p, nrow=n, byrow=T)
uXb = Xb / rowSums(Xb^2)
mean.mat.b = matrix(0, ncol=p,nrow=nboot)
for(i in 1:nboot){
  mean.mat.b[i,] = apply(uXb[sample(1:n,n,replace=T),], 2, mean)
}
dep.u = mdepth.MhD(mean.mat.b, mean.mat.b)$dep
dep.ub = mdepth.MhD(t(as.matrix(uX0bar)), mean.mat.b)$dep
(pval.u = sum(dep.ub > dep.u)/nboot)
