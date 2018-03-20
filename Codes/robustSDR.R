setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(rrcov)
library(fda.usc)
library(ldr)

n = 80
V = matrix(0, ncol=2, nrow=n)
for(i in 1:n){
  ri = sample(1:4, 1)
  if(ri==1){
    V[i,] = c(-1,runif(1,-1,1))
  }
  else if(ri==2){
    V[i,] = c(runif(1,-1,1), 1)
  }
  else if(ri==3){
    V[i,] = c(1,runif(1,-1,1))
  }
  else{
    V[i,] = c(runif(1,-1,1), -1)
  }
}

plot(V, pch=19, cex=.3)

p = 1e2
d = 2
X = V %*% matrix(rnorm(d*p), nrow=d) + matrix(rnorm(n*p), ncol=p)

X1 = X
X1[1:5,1] = X1[1:5,1] + 1e5

par(mfrow=c(1,2))
z = PcaLocantore(t(X1), k=2)
plot(z@loadings, cex=.5, pch=19)

z = PcaRank(t(X1), k=2)
plot(z@loadings, cex=.5, pch=19)
par(mfrow=c(1,1))

# solver function:
rsir = function(X,y,d,nslices){
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)

  sy = ldr.slices(y, nslices=nslices)$slice.indicator
  uy = unique(sy)
  luy = length(uy)
  
  # get unique mean levels
  u.mu = matrix(0, ncol=p, nrow=as.numeric(luy))
  for(i in 1:length(uy)){
#    u.mu[i,] = apply(X[which(sy == uy[i]),], 2, median)
    u.mu[i,] = spatial.median(X[which(sy == uy[i]),], delta=1e-5)$mu
    
  }
  Xbar = spatial.median(X, delta=1e-5)$mu
  u.mu = u.mu - matrix(Xbar, ncol=p, nrow=luy, byrow=T)
  
  # get Gamma
  pcamod = PcaRank(X)
  Gamma = pcamod@loadings[,1:d]
  
  # Get signal part of X
  muhat = X
  for(i in 1:n){
    muhat[i,] = u.mu[which(uy == sy[i]),]
  }
  muhat = muhat %*% Gamma %*% t(Gamma) +
    matrix(Xbar, ncol=p, nrow=n, byrow=T)
  sigma2hat = mean((apply(X - muhat, 2, mad))^2)

  a = list(X=X, y=y, d=d, Gammahat=Gamma, muhat=muhat, sigma2hat=sigma2hat)
  class(a) = "rsirmod"
  a
}

sir = function(X,y,d,nslices){
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  
  sy = ldr.slices(y, nslices=nslices)$slice.indicator
  uy = unique(sy)
  luy = length(uy)
  
  # get unique mean levels
  u.mu = matrix(0, ncol=p, nrow=luy)
  for(i in 1:length(uy)){
    u.mu[i,] = apply(X[which(sy == uy[i]),], 2, mean)
  }
  u.mu = u.mu - matrix(apply(X,2,mean), ncol=p, nrow=luy, byrow=T)
  
  # get Gamma
  pcamod = PcaClassic(X)
  Gamma = pcamod@loadings[,1:d]
  
  # Get signal part of X
  muhat = X
  for(i in 1:n){
    muhat[i,] = u.mu[which(uy == sy[i]),]
  }
  muhat = muhat %*% Gamma %*% t(Gamma) +
    matrix(apply(X,2,mean), ncol=p, nrow=n, byrow=T)
  sigma2hat = mean((apply(X - muhat, 2, sd))^2)
  
  a = list(X=X, y=y, d=d, Gammahat=Gamma, muhat=muhat, sigma2hat=sigma2hat)
  class(a) = "rsirmod"
  a
}

# predict
predict.rsir = function(rsirmod, newx){
  R.X = with(rsirmod, X %*% Gammahat)
  R.newx = with(rsirmod, newx %*% Gammahat)
  n = nrow(rsirmod$X)
  
  pred = rep(0, nrow(newx))
  for(i in 1:nrow(newx)){
    R.idiff = matrix(as.numeric(R.newx[i,]), ncol=rsirmod$d, nrow=n, byrow=T) - R.X
    wi = exp(- rowSums(R.idiff^2)/(2*rsirmod$sigma2hat))
    pred[i] = sum(rsirmod$y*wi) / sum(wi)
  }
  pred
}

# repeat for 100 datasets
p.vec = c(5,10,25,50,75,100,125,150)
n = 200

err.mat = matrix(0, ncol=2, nrow=length(p.vec))
pb = txtProgressBar(0,length(p.vec))
for(i in 1:length(p.vec)){
  p = p.vec[i]
  err = matrix(0,nrow=1e2,ncol=2)
  for(nsim in 1:1e2){
    # generate data
    y = rnorm(2*n)
    vy = y + y^2 + y^3
    Gamma = rep(1,p)
    X = outer(vy, Gamma) + 5*matrix(rnorm(2*n*p), ncol=p)
    X[1:10,1:(p/5)] = X[1:10,1:(p/5)]+100
    
    # fit models
    mod = sir(X[1:n,],y[1:n],d=1,nslices=20)
    rmod = rsir(X[1:n,],y[1:n],d=1,nslices=20)
    
    # predict for new data
    pred.vec = predict.rsir(mod, X[(n+1):(2*n),])
    err[nsim,1] = mean((pred.vec-y[(n+1):(2*n)])^2)
    pred.vec = predict.rsir(rmod, X[(n+1):(2*n),])
    err[nsim,2] = mean((pred.vec-y[(n+1):(2*n)])^2)
  }
  
  err.mat[i,1] = mean(err[,1], na.rm=T)
  err.mat[i,2] = mean(err[,2], na.rm=T)
  setTxtProgressBar(pb, i)
}
close(pb)

err.mat
pdf('SDRcomparison_noout.pdf', width=5, height=5)
plot(p.vec, err.mat[,2], lty=1, lwd=2, type='b', xlab="p", ylab="Prediction errors", ylim=c(0,1))
lines(p.vec, err.mat[,1], lty=2, lwd=2, type='b')
legend('topright', c("SDR","robust SDR"), lty=c(2,1), lwd=2)
dev.off()

## PFC simulation: application on big mac data
data(bigmac)
bigmacr = as.matrix(scale(bigmac[,-1]))

ipred = rep(0 ,nrow(bigmac))
for(i in 1:nrow(bigmac)){
  imod = rsir(X=bigmacr[-i,], y=bigmac[-i,1], d=1, nslices=5)
  
  ipred[i] = predict.rsir(imod, t(as.matrix(bigmacr[i,])))
}

mean((ipred - bigmac[,1])^2)
  
