library(plot3D)
library(fda.usc)
library(scales)

X = matrix(rnorm(1e3), ncol=2) %*% sqrt(diag(c(1,1/2)))

# make grid of points
pts = seq(-3, 3, by=.2)
lengrid = length(pts)
xcoord = rep(pts, rep(lengrid,lengrid))
ycoord = rep(pts, lengrid)
xygrid = cbind(xcoord,ycoord)
rm(xcoord,ycoord)

pdf('depthplot.pdf', width=6, height=5)
p = persp3D(pts, pts, matrix(mdepth.RP(xygrid,X)$dep, nrow=lengrid, byrow=T),
      xlab="x1", ylab="x2", zlab="Depth", zlim=c(-1,1),
      col=alpha("lightblue", .2),expand = 0.5,shade = 0.2,
      theta=30, phi=10, border=alpha("black", .2), contour=TRUE)

points3D(X[,1], X[,2], rep(-1,nrow(X)),
         pch=19, cex=.5, pch=.3, col='red', add=TRUE)
dev.off()