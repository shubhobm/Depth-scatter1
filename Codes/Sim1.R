## PQDSim: Initial simulations of Projection quantile depth and comparisions with projection depth
## v4: Multivariate rank-based PCA simulations

## Functions
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

## function to calculate weighted projection quantile depth
EPQD1 = function(X, grid, nu=1e3){
  
  p = ncol(X)
  b = apply(X, 2, median)
  X0 = X - ones(nrow(X),1) %*% b
  grid0 = grid - ones(nrow(grid),1) %*% b
  
  ## get matrix of weighted PQDs for all points
  npt = dim(grid)[1]
  Fuxu.mat = matrix(0, nrow=npt, ncol=nu)
  
  # loop over nu pts on unit circle then take max
  for(iu in 1:nu){
    u = as.matrix(rnorm(p)); u = u/sqrt(sum(u^2))
    uecdf = ecdf(X0%*%u)
    Fuxu.mat[,iu] = uecdf(grid0%*%u)
  }
  EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  
  return(cbind(grid,EPQD.vec))
  
}

wEPQD1 = function(X, grid, sig, nu=1e3){
  
  p = ncol(X)
  b = apply(X, 2, median)
  X0 = X - ones(nrow(X),1) %*% b
  grid0 = grid - ones(nrow(grid),1) %*% b
  
  ## get matrix of weighted PQDs for all points
  npt = dim(grid)[1]
  Fuxu.mat = matrix(0, nrow=npt, ncol=nu)
  
  # loop over nu pts on unit circle then take max
  for(iu in 1:nu){
    u = as.matrix(rnorm(p)); u = u/sqrt(sum(u^2))
    I.minus.Pu = diag(p) - u%*%t(u)
    
    Xuperp = X0 %*% I.minus.Pu
    scaled.perp = sqrt(Xuperp^2 %*% ones(ncol(X),1))
    #w = ifelse(scaled.perp>sig, 0, 1)
    #w = sig*exp(-scaled.perp/sig)
    w = dnorm(scaled.perp, sd=sig)
    #w = dcauchy(Xuperp, scale=sig)
    uecdf = ecdf(w * (X0%*%u))
    
    gridperp = grid0 %*% I.minus.Pu
    scaled.gridperp = sqrt(gridperp^2 %*% ones(ncol(X),1))
    #wu = ifelse(scaled.gridperp>sig, 0, 1)
    #wu = sig*exp(-scaled.gridperp/sig)
    wu = dnorm(scaled.gridperp, sd=sig)
    #wu = dcauchy(sqrt(apply(xygrid^2,1,sum) - xygrid.u^2), scale=sig)
    Fuxu.mat[,iu] = uecdf(wu * (grid0%*%u))
  }
  EPQD.vec = 1/(1+apply(abs(Fuxu.mat-.5), 1, max))
  
  return(cbind(grid,EPQD.vec))
  
}

pcarank = function(X, ...){
  X = as.matrix(X)
  d1 = EPQD1(X, X)
  X1 = X
  for(i in 1:nrow(X1)){
    X1[i,] = X1[i,]/sqrt(sum(X1[i,]^2))
  }
  Xrank = X1 * (1/d1[,3]-1)
  princomp(Xrank, ...)
}

pcarank1D = function(X, ...){
  X1 = apply(X,2,rank)
  princomp(X1, ...)
}

## Empirically calculate PQD with grid search
# Bivariate normal mixture
n = 1e3
set.seed(120214)
# x = rnorm(n)
# X = cbind(x, 2*x)
# X = scale(X, scale=F)
Gamma = matrix(c(1,-1,1,1), nrow=2)/sqrt(2)
sig = Gamma %*% diag(c(25,1)) %*% t(Gamma)
X = my.mvrnorm(n, mu=c(0,0), Sig=2*sig)

# x = runif(n)
# X = cbind(x, x+rnorm(n, sd=.001))
# X = scale(X, scale=F)

# PCA before contamination
(p.pure <- princomp(X)); p.pure$loadings

d1 = EPQD1(X, X)
X1 = X
for(i in 1:nrow(X1)){
  X1[i,] = X1[i,]/sqrt(sum(X1[i,]^2))
}
Xrank = X1
Xrank = X1 * (max(d1[,3]) - d1[,3])
(p.rank = princomp(Xrank)); p.rank$loadings

par(mfrow=c(1,2))
plot(X, pch=19, cex=.5)
plot(Xrank, pch=19, cex=.5)
par(mfrow=c(1,1))

# after contamination
pure = 600
X[(pure+1):n,] = X[(pure+1):n,] +
  matrix(sample(c(10:120, -120:-10), (n-pure)*2, replace=T), ncol=2)
(p.impure <- princomp(X)); p.impure$loadings

d1 = EPQD1(X, X)
X1 = X
for(i in 1:nrow(X1)){
  X1[i,] = X1[i,]/sqrt(sum(X1[i,]^2))
}
Xrank = X1 * (1/d1[,3]-1)
(p.rank = princomp(Xrank)); p.rank$loadings

## data example
data("heptathlon", package = "HSAUR")
hep1 = as.matrix(heptathlon[,-8])
(ph = princomp(hep1, center=T, scale=F)); ph$loadings

(ph.rank = pcarank(hep1, scale=F)); ph.rank$loadings

(ph.rank1D = pcarank1D(hep1, scale=F)); ph.rank1D$loadings

require(rrcov)
ph1 = PcaHubert(heptathlon[,-8])
ph1@loadings

# compare with robust pca
require(rrcov)
p2 = PcaHubert(hbk)
p2@loadings

p2.rank = pcarank(hbk)
p2.rank$loadings
p2.rank$sdev

hbk1 = hbk[-(1:14),-4]
p21 = PcaHubert(hbk1)
p21@loadings

p21.rank = pcarank(hbk1)
p21.rank$loadings
p21.rank$sdev

p21.pc = princomp(hbk1)
p21.pc$loadings
