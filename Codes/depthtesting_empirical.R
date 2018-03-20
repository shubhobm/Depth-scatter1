## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')

library(fda.usc)

get.power = function(n,p,nrep,mu,cov,rho=0.5){
  
  ## construct output matrix
  out.mat = matrix(0, nrow=length(n)*length(p), ncol=6)
  out.mat[,1] = rep(p,rep(length(n),length(p)))
  out.mat[,2] = rep(n, length(p))
  
  ## get rows of output matrix, each by nrep simulations
  loopfun = function(j){
    is.reject = rep(0,4)
    
    # generate data
    X1 = my.mvrnorm(ni, mu=rep(0,pi), Sig=Sigma)
    #X2 = my.mvrnorm(ni, mu=rep(mu,pi), Sig=9*Sigma)
    X2 = matrix(rcauchy(pi*ni), ncol=pi)
    ind = rbinom(ni, 1, 0.5)
    X = mu + ind*X1 + (1-ind)*X2
    
    # Hotelling Tsq
    Xbar = apply(X, 2, mean)
    T2 = ni * t(Xbar) %*% solve(cov(X)) %*% Xbar
    is.reject[1] = (1 - pchisq(T2, pi) < 0.05)
    
    ## Sign-based test
    uX = X / sqrt(rowSums(X^2))
    uXbar = apply(uX, 2, mean)
    uT2 = ni * t(uXbar) %*% solve(cov(uX)) %*% uXbar
    is.reject[2] = (sum(uDist > as.numeric(uT2))/1e4 < 0.05)
    
    ## Depth-rank based test
    dX = mdepth.RP(X,X)$dep
#   Xc = X - apply(X,2,median)
#    uXc = Xc / sqrt(rowSums(Xc^2))
    rX = uX * (max(dX)-dX)
    rXbar = apply(rX, 2, mean)
    rT2 = ni * t(rXbar) %*% solve(cov(rX)) %*% rXbar
    is.reject[3] = (sum(rDist > as.numeric(rT2))/1e4 < 0.05)

    ## Depth-based test
    Xd = uX * dX
    dXbar = apply(Xd, 2, mean)
    dT2 = ni * t(dXbar) %*% solve(cov(Xd)) %*% dXbar
    is.reject[4] = (sum(dDist > as.numeric(dT2))/1e4 < 0.05)
    #c(T2,uT2,rT2)
    
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
    
    ## get null distributions
    #DSigma = diag(sqrt(eigen(Sigma)$values))
    #Samp = matrix(rnorm(1e4*pi), ncol=pi) %*% DSigma
    Samp1 = my.mvrnorm(1e4, mu=rep(0,pi), Sig=Sigma)
    Samp2 = matrix(rcauchy(pi*1e4), ncol=pi)
    ind = rbinom(1e4, 1, 0.5)
    Samp = ind*Samp1 + (1-ind)*Samp2
    uSamp = Samp / sqrt(rowSums(Samp^2))
    uDist = diag(uSamp %*% solve(cov(uSamp)) %*% t(uSamp))
    
    depSamp = mdepth.RP(Samp,Samp)$dep
    rSamp = uSamp * (max(depSamp) - depSamp)
    rDist = diag(rSamp %*% solve(cov(rSamp)) %*% t(rSamp))
    
    dSamp = uSamp * depSamp
    dDist = diag(dSamp %*% solve(cov(dSamp)) %*% t(dSamp))
        
    power.vec = lapply(1:nrep, loopfun)
    #power.vec = mclapply(1:nrep, loopfun, mc.cores=4) # parallel version of lapply
    power.vec = apply(matrix(unlist(power.vec), ncol=4, byrow=T), 2, mean)
    out.mat[irep,3:6] = power.vec
  }
  
  out.frame = data.frame(out.mat)
  names(out.frame) = c("p","n", "T2","Sign","Rank","Depth")
  out.frame
}

n = c(50,100,200,500)
p = c(5,10,20)
nrep = 1e2

get.power(n, p, nrep, mu=0.2, cov="AR1", rho=0.8)
get.power(n, p, nrep, mu=0, cov="compound", rho=0.2)
#get.power(n, p, nrep, mu=0.2, cov="srivatsa")

