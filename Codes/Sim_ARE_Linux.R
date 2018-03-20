## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

#setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(parallel)

## Functions

## simulate for normal dist
FSE.norm = function(n, p, iter=1e3, ncores=detectCores()){
  set.seed(12182014)
  v = c(rep(0,p-1), 1)
  lam = 1:p
  Sigma = diag(lam)
  
  MSE.mat = matrix(0, nrow=length(n), ncol=8)
  for(i in 1:length(n)){
    
    # function to compute stuff 1000 times for a given n
    loopfun = function(j){
      source('misc_functions.R')
      require(fastM)
      require(fda.usc)
      iv = rep(0,9)
      
      # get sample and construct sign matrix
      iX = matrix(rnorm(p*n[i]), ncol=p) %*% sqrt(Sigma)
      iXnorm = sqrt(iX^2 %*% rep(1,p))
      iS = iX / (iXnorm %*% rep(1,p))
      
      # PCA on original sample
      iP = princomp(iX)
      iv[1] = abs(sum(v * iP$loadings[,1]))
      
      # PCA on SCM
      iPsign = princomp(iS)
      iv[2] = abs(sum(v * iPsign$loadings[,1]))
      
      # PCA on Tyler's cov matrix
      T = TYLERshape(iX)$Sigma
      iv[3] = abs(sum(v * eigen(T)$vectors[,1]))
      
      # PCA on DCM and depth-weighted tyler's scatter
      # Tukey's depth
      idep = mdepth.HS(iX,iX)$dep
      #idep = EPQD(iX,iX)[,p+1]
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[4] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[5] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      # Mahalanobis depth
      idep = mdepth.MhD(iX,iX)$dep
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[6] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[7] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      # Projection depth
      idep = mdepth.RP(iX,iX)$dep
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[8] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[9] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      iv
    }
    
    # parallel code: compute MSE elements iter times
    #system.time(eff.v <- foreach(j=1:iter) %dopar% loopfun(j))
    system.time(eff.v <- mclapply(1:iter, loopfun, mc.cores=ncores))
    
    # get MSE and return
    eff.v = matrix(unlist(eff.v), ncol=9, byrow=T)
    (MSE.vec = apply(eff.v, 2, function(x) mean(acos(x)^2)))
    MSE.mat[i,] = MSE.vec[1]/MSE.vec[-1]
  }
  
  MSE.mat
}

## simulate for t-distn
FSE.t = function(n, p, df, iter=1e3, ncores=detectCores()){
  set.seed(12182014)
  v = c(rep(0,p-1), 1)
  lam = 1:p
  Sigma = diag(lam)
  
  MSE.mat = matrix(0, nrow=length(n), ncol=8)
  for(i in 1:length(n)){
    
    # function to compute stuff 1000 times for a given n
    loopfun = function(j){
      source('misc_functions.R')
      require(fastM)
      require(fda.usc)
      require(mvtnorm)
      iv = rep(0,9)
      
      # get sample and construct sign matrix
      iX = rmvt(n[i], sigma=Sigma, df=df)
      iXnorm = sqrt(iX^2 %*% rep(1,p))
      iS = iX / (iXnorm %*% rep(1,p))
      
      # PCA on original sample
      iP = princomp(iX)
      iv[1] = abs(sum(v * iP$loadings[,1]))
      
      # PCA on SCM
      iPsign = princomp(iS)
      iv[2] = abs(sum(v * iPsign$loadings[,1]))
      
      # PCA on Tyler's cov matrix
      T = TYLERshape(iX)$Sigma
      iv[3] = abs(sum(v * eigen(T)$vectors[,1]))
      
      # PCA on DCM and depth-weighted tyler's scatter
      # Tukey's depth
      idep = mdepth.HS(iX,iX)$dep
      #idep = EPQD(iX,iX)[,p+1]
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[4] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[5] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      # Mahalanobis depth
      idep = mdepth.MhD(iX,iX)$dep
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[6] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[7] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      # Projection depth
      idep = mdepth.RP(iX,iX)$dep
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[8] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[9] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      iv
    }
    
    # parallel code: compute MSE elements iter times
    #system.time(eff.v <- foreach(j=1:iter) %dopar% loopfun(j))
    system.time(eff.v <- mclapply(1:iter, loopfun, mc.cores=ncores))
    
    # get MSE and return
    eff.v = matrix(unlist(eff.v), ncol=9, byrow=T)
    (MSE.vec = apply(eff.v, 2, function(x) mean(acos(x)^2)))
    MSE.mat[i,] = MSE.vec[1]/MSE.vec[-1]
  }
  
  MSE.mat
}


# n.vec = c(20,40,60,100,150,200,250,300)
# system.time(norm.table <- FSE.norm(n.vec, 2, 1e3))
# m5 = 1
# plot(norm.table[,1]~n.vec, type='b', lty=2, ylim=c(0,m5+.1),
#      xlab="sample size", ylab="Efficiency")
# lines(norm.table[,2]~n.vec, type='b', lty=1)
# lines(norm.table[,3]~n.vec, type='b', lty=2, pch=2, col='red')
# lines(norm.table[,4]~n.vec, type='b', lty=1, pch=2, col='red')
# lines(norm.table[,5]~n.vec, type='b', lty=2, pch=3, col='blue')
# lines(norm.table[,6]~n.vec, type='b', lty=1, pch=3, col='blue')
# lines(norm.table[,7]~n.vec, type='b', lty=2, pch=4, col='darkgreen')
# lines(norm.table[,8]~n.vec, type='b', lty=1, pch=4, col='darkgreen')


# Large sample efficiencies
system.time(ARE.t5 <- FSE.t(1e3, 2, df=5, 1e4, ncores=16))
system.time(ARE.t6 <- FSE.t(1e3, 2, df=6, 1e4, ncores=16))
system.time(ARE.t10 <- FSE.t(1e3, 2, df=10, 1e4, ncores=16))
system.time(ARE.t15 <- FSE.t(1e3, 2, df=15, 1e4, ncores=16))
system.time(ARE.t25 <- FSE.t(1e3, 2, df=25, 1e4, ncores=16))
system.time(ARE.norm <- FSE.norm(1e3, 2, 1e4, ncores=16))
ARE.table = rbind(ARE.t5, ARE.t6, ARE.t10, ARE.t15, ARE.t25, ARE.norm)
ARE.table  =ARE.table[,c(1,2,3,5,7,4,6,8)]
write.csv(format(round(ARE.table, 2), nsmall=2), "tableARE.txt")
