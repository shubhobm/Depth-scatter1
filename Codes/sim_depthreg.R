rm(list=ls())

## Function to generate multivariate normal observations
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

## simulation setup
n = 1e3; p = 4
v.star = 1
beta.star = 1:4
Sigma = diag(p)

X = my.mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma)
y = my.mvrnorm(n=1, mu=rep(0,n)+X%*%beta.star, Sigma=v.star*diag(n))

npure = 600
y[(npure+1):n] = y[(npure+1):n] + sample(60:109, (n-npure), replace=T)

beta.mat = matrix(0, nrow=n, ncol=p+1)
idep = rep(1,n)

iter=1; maxit=100
for(i in 1:n){
  ilm = lm(t(y)~X, weights=idep)
  beta.mat[i,] = ilm$coef
  ires = residuals(ilm)
  iPhi = pnorm(ires)
  idep = iPhi*(1-iPhi)
}

par(mfrow=c(2,2))
for(i in 2:5){
  plot(1:n, beta.mat[,i], type="l")
}
par(mfrow=c(1,1))

col.vec = c(rep("green", npure), rep("red", n-npure))
plot(idep, col=col.vec, pch=19, cex=.5)

plot(ilm$res~ilm$fit, col=col.vec, pch=19)
