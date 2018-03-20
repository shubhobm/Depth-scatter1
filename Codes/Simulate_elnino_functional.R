setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')
library(rrcov)
library(fda.usc)
library(splines)
library(pracma)
library(robustbase)

npfda.elnino <- read.table(
  "C:/Study/My projects/Depth-scatter/Data/npfda-elnino.dat", quote="\"")
mydata = npfda.elnino
# mydata = Octane.X
n = nrow(mydata)
p = ncol(mydata)

# take a B-spline basis and orthonormalize it
B = bs(1:p, knots=seq(0,p-1,by=2), degree=3)
B1 = gramSchmidt(B)$Q

# transform data... diff between two bins is 1
Xt = as.matrix(mydata[,-1]) %*% B1[-1,]
muhat = spatial.median(Xt, delta=1e-5)$mu

set.seed(08162016)
source('misc_functions.R')
ncomp = 2

nrep = 1e3
ndata = 1e2
out.mat = matrix(0, nrep, 3)

## iterate nrep times
pb = txtProgressBar(0, nrep)
for(rep in 1:nrep){
  pcarank = PcaRank(Xt, k=ncomp, proj=1e2)
  mu = muhat + colMeans(pcarank@scores %*% t(pcarank@loadings))

  # generate random observations
  V = c(1,1,1,1,10)
  nbasis = ncol(B)

  ## add noise to loadings
  B.mat = matrix(0, nrow=ndata, ncol=nbasis)
  for(i in 1:ndata){
    B.mat.j = matrix(0, nrow=nbasis, ncol=ncomp)
    for(j in 1:ncomp){
      B.mat.j[,j] = my.mvrnorm(n=1, mu=pcarank@loadings[,j], Sigma=diag(rep(V[j],nbasis)))
    }
    B.mat[i,] = muhat + colMeans(pcarank@scores %*% t(B.mat.j))
  }
  
  ## add white noise to each element
  B.mat = B.mat + matrix(rnorm(ndata*nbasis, 0, 1e-1), ndata, nbasis)
  
  ## add comtamination
  B.mat[2:6,2:6] = B.mat[2:6,2:6] + 1e2
  
  out.mat[rep,] = c(sqrt(mean((mu - spatial.median(B.mat, delta=1e-5)$mu)^2)),
                  sqrt(mean((mu - depth.median(B.mat, delta=1e-5)$mu)^2)),
                  sqrt(mean((mu - colMeans(B.mat))^2)))
  setTxtProgressBar(pb,rep)
}
close(pb)


Xhat = B.mat %*% t(B1)
# Xhat[,1] = as.numeric(mydata[,1])
plot(as.numeric(Xhat[1,-1]), col=grey(.6), type='l',xaxt='n',
     ylab="Sea surface temperature", xlab="Month",
     ylim=(c(min(Xhat), max(Xhat+.1))))
axis(side=1, at=1:12, labels=c('Jun','Jul','Aug','Sep','Oct','Nov','Dec',"Jan","Feb",'Mar','Apr','May'))
for(i in 2:n){
  lines(as.numeric(Xhat[i,-1]), col=grey(.6))
}
lines((B1 %*% mu)[-1], lwd=2)
lines((B1 %*% spatial.median(B.mat, delta=1e-5)$mu)[-1], lwd=2, col=adjustcolor("red", alpha.f=.5))
lines((B1 %*% depth.median(B.mat, delta=1e-5)$mu)[-1], lwd=2, col=adjustcolor("blue", alpha.f=.5))
lines((B1 %*% colMeans(B.mat))[-1], lwd=2, col=adjustcolor("green", alpha.f=.5))


apply(out.mat, 2, mean) # mean error for spatial median, depth median and mean
apply(out.mat, 2, sd) # sd error for spatial median, depth median and mean
