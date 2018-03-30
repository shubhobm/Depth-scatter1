## Sim_influenceplot: Plots to compare influence functions for DCM and SCM

setwd("d:/Study/My projects/Depth-scatter1/Codes")
rm(list=ls());
source('misc_functions.R')

## Functions
# DCM plots for different types of depth
plot.IFnorm = function(X, grid, depth, ...){
  require(fda.usc)
  
  # get depth values
  if(depth=='HS'){
    DX = mdepth.HS(X,X)$dep
  }
  else if(depth=='MhD'){
    DX = mdepth.MhD(X,X)$dep
  }
  else {
    DX = mdepth.RP(X,X)$dep
  }
    
  # get eigenvalues of DCM
  DX = max(DX) - DX
  lamDS = colMeans((DX^2/sum.lam.Zsq) * lam.Zsq)
  
  # get htped at grid points
  if(depth=='HS'){
    Dgrid = mdepth.HS(grid,X)$dep
  }
  else if(depth=='MhD'){
    Dgrid = mdepth.MhD(grid,X)$dep
  }
  else {
    Dgrid = mdepth.RP(grid,X)$dep
  }
  Dgrid = max(Dgrid) - Dgrid
  
  # get norms of influence fns for eigenvectors
  mult = sqrt(lam[1]*lam[2])/(lamDS[1] - lamDS[2])
  IFnorm.D = abs(mult * xygrid[,1] * xygrid[,2] * Dgrid^2 / diag(xygrid %*% Sig %*% t(xygrid)))
  
  # plot result
  persp(pts, pts, matrix(IFnorm.D, nrow=lengrid, byrow=T), ...)
}

## Scatterplot of data and D-rank
require(fda.usc)
n = 1e3
set.seed(120214)
Gamma = matrix(c(1,-1,1,1), nrow=2)/sqrt(2)
sig = Gamma %*% diag(c(9,1)) %*% t(Gamma)
X = my.mvrnorm(n, mu=c(0,0), Sig=2*sig)
X = rbind(X, matrix(c(-5,-9, -7,-9, 9,9), ncol=2, byrow=T))
uX = X / sqrt(rowSums(X^2))
dX = mdepth.HS(X,X)$dep
Xrank = uX * (max(dX) - dX)

pdf('signs.pdf', width=7, height=4)
par(mfrow=c(1,2))
plot(X, pch=19, cex=.5)
plot(uX, pch=19, cex=.5)
par(mfrow=c(1,1))
dev.off()

pdf('ranks.pdf', width=7, height=4)
par(mfrow=c(1,2))
plot(X[1:1000,], pch=19, cex=.5)
points(X[1001:1003,], col='red', pch=19, cex=1)
plot(Xrank[1:1000,], pch=19, cex=.5)
points(Xrank[1001:1003,], col='red', pch=19, cex=1)
par(mfrow=c(1,1))
dev.off()

pdf('signs_and_ranks.pdf', width=8, height=3)
par(mfrow=c(1,3))
plot(X[1:1000,], pch=19, cex=.5)
points(X[1001:1003,], col='red', pch=19, cex=1)
plot(uX[1:1000,], pch=19, cex=.5)
points(uX[1001:1003,], col='red', pch=19, cex=1)
plot(Xrank[1:1000,], pch=19, cex=.5)
points(Xrank[1001:1003,], col='red', pch=19, cex=1)
par(mfrow=c(1,1))
dev.off()

## Influence function plots
lam = c(2,1)
Sig = diag(lam)
Z = matrix(rnorm(1e3),ncol=2)
X = Z %*% sqrt(Sig)

# make grid of points
pts = seq(-3, 3, by=.2)
lengrid = length(pts)
xcoord = rep(pts, rep(lengrid,lengrid))
ycoord = rep(pts, lengrid)
xygrid = cbind(xcoord,ycoord)
rm(xcoord,ycoord)

default = par()
par(mfrow=c(3,2), mai=rep(.5,4))
# Sample covariance matrix
r = sqrt(rowSums(xygrid^2))
Ugrid = xygrid / r
IFnorm = abs(Ugrid[,1] * Ugrid[,2] * sqrt(lam[1]*lam[2])/(lam[1] - lam[2]))

IFnorm.Sig = r^2*IFnorm
persp(pts, pts, matrix(IFnorm.Sig, nrow=lengrid, byrow=T),
      main="(a)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

# Influence fn plot for SCM
# get eigenvalues of SCM
lam.Zsq = (Z * Z) %*% Sig
sum.lam.Zsq = rowSums(lam.Zsq)
lamS = colMeans(lam.Zsq / sum.lam.Zsq)

# get norms of influence fns for eigenvectors
multS = sqrt(lam[1]*lam[2])/(lamS[1] - lamS[2])
IFnorm.S = sqrt(abs(multS * Ugrid[,1] * Ugrid[,2]))

# plot result
persp(pts, pts, matrix(IFnorm.S, nrow=lengrid, byrow=T),
      main="(b)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

# Influence fn plot for Tyler's scatter matrix
IFnorm.tyler = 4*IFnorm
persp(pts, pts, matrix(IFnorm.tyler, nrow=lengrid, byrow=T),
      main="(c)", xlab="x1", ylab="x2", zlab="IF(x0)",
      ticktype="detailed", nticks=3,
      theta=45, phi=45, col=gray(.9), border=gray(.3))

# Influence fn plot for DCM

# Tukey's halfspace depth
plot.IFnorm(X, xygrid, 'HS',
            main="(d)", xlab="x1", ylab="x2", zlab="IF(x0)",
            ticktype="detailed", nticks=3,
            theta=45, phi=45, col=gray(.9), border=gray(.3))

# Mahalanobis depth
plot.IFnorm(X, xygrid, 'MhD',
            main="(e)", xlab="x1", ylab="x2", zlab="IF(x0)",
            ticktype="detailed", nticks=3,
            theta=45, phi=45, col=gray(.9), border=gray(.3))

# Zuo's projection depth
plot.IFnorm(X, xygrid, 'RP',
            main="(f)", xlab="x1", ylab="x2", zlab="IF(x0)",
            ticktype="detailed", nticks=3,
            theta=45, phi=45, col=gray(.9), border=gray(.3))

par(default)

# persp(pts, pts, matrix(mdepth.HS(xygrid,X)$dep, nrow=lengrid, byrow=T),
#       main="(a)", xlab="x1", ylab="x2", zlab="IF(x0)",
#       ticktype="detailed", nticks=3,
#       theta=45, phi=45, col=gray(.9), border=gray(.3))

