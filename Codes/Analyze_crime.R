## Comm_analyze: analysis of communities data
rm(list=ls())
setwd("C:/Study/My projects/Depth-scatter/Codes")
source('misc_functions.R')

##########
# Load packages
##########
library(fda.usc)
library(rrcov)

##########
# Load data
##########
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

##########
# Implementation
##########

## wil do vanilla and rank PCA, then
## principal components regression on increasing number of PCs

## vanilla PCA
pcamod = PcaClassic(commdata.X)
scores = pcamod@scores

# rank PCA
pca.rank = PcaRank(commdata.X)
scores.rank = pca.rank@scores

par(mfrow=c(1,2))
barplot(pcamod@eigenvalues[1:10]/sum(pcamod@eigenvalues),
        ylim=c(0,.4), main="Normal PCA")
barplot(pca.rank@eigenvalues[1:10]/sum(pca.rank@eigenvalues),
        ylim=c(0,.4), main="(Projection-) Depth PCA")
par(mfrow=c(1,1))

n = nrow(commdata.X)
p = ncol(commdata.X)
rsq.mat = matrix(0, nrow=p, ncol=2)
for(npc in 1:p){
  # PC regression on vanilla PCA scores
  lm.PC = lm(1/Y^3~scores[,1:npc]) # inverse cube transformation chosen by boxcox
  
  # PC-regression on rank PCA scores
  lm.PCrank = lm(1/Y^3~scores.rank[,1:npc])
  
  rsq.mat[npc,] = c(summary(lm.PC)$r.squared, summary(lm.PCrank)$r.squared)
}

plot(1:p, rsq.mat[,1], type='l', ylim=c(0.2,.8), lwd=2,
     xlab="no. of PCs", ylab="R^2")
lines(1:p, rsq.mat[,2], type='l', lty=2, lwd=2)
legend("bottomright", c("with normal PCs", "with depth-PCs"), lty=c(1,2), lwd=2)
