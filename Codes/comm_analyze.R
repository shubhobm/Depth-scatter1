## Comm_analyze: analysis of communities data
rm(list=ls())
setwd("C:/Study/My projects/Depth-scatter/Codes")
commdata = read.csv("../Data/communities.data.txt", header=F)
NAcols = c(102:118, 122:125, 127) # columns with most NA entries
commdata.X = commdata[,-c(1:5, NAcols)]

# convert all columns into numeric
for(i in 1:ncol(commdata.X)){
  if(class(commdata.X[,i]) != "numeric"){
    commdata.X[,i] = as.numeric(paste(commdata.X[,i]))
  }
}

Y = commdata.X[,101]+1
commdata.X = data.frame(scale(commdata.X[,-101]))
Y = Y[complete.cases(commdata.X)]
commdata.X = commdata.X[complete.cases(commdata.X),]

## vanilla PCA
pcamod = princomp(commdata.X)
plot(pcamod)
scores = pcamod$scores
pairs(scores[,1:3], pch=19, cex=.5)

# PC-regression
X.PC = as.matrix(commdata.X) %*% pcamod$loadings[,1:5]
# best transformation: boxcox
boxcox(Y~X.PC, seq(-5,5,.1))
lm.PC = lm(1/Y^3~X.PC)
summary(lm.PC)

par(mfrow=c(2,2))
plot(lm.PC)
par(mfrow=c(1,1))

## rank PCA
n = nrow(commdata.X)
p = ncol(commdata.X)
norms = sqrt(commdata.X^2 %*% rep(1,p))
signs = commdata.X / (norms %*% rep(1,p))

# calculate depth
require(fda.usc)
depths = mdepth.RP(commdata.X, commdata.X)$dep
depths = max(depths) - depths
commdata.rank = signs * depths

pca.rank = princomp(commdata.rank)
plot(pca.rank)

scores.rank = pca.rank$scores
pairs(scores.rank[,1:3], pch=19, cex=.5)

# PC-regression
X.PCrank = as.matrix(commdata.rank) %*% pca.rank$loadings[,1:5]

# best transformation: boxcox
boxcox(Y~X.PCrank, seq(-5,5,.1))

lm.PCrank = lm(1/Y^3~X.scores.rank[,1:5])
summary(lm.PCrank)

par(mfrow=c(2,2))
plot(lm.PCrank)
par(mfrow=c(1,1))