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

# rain = read.csv("C:/Study/My projects/Climate-indian monsoon/data/rain_1973_2013_test.csv")
# rainsmall = aggregate(PRCP ~ year+STATION_NAME, data=rain, FUN=median)
# rainwide = reshape(rainsmall, timevar="year", idvar="STATION_NAME", direction="wide")

# colvec = rep(grey(.5), n)
# colvec[c(6,25,26,36:39)] = "black"
# plot(as.numeric(Octane.X[1,]), type='l', col=colvec[1])
# for(i in 2:n){
#   lines(as.numeric(Octane.X[i,]), col=colvec[i])
# }

mydata = npfda.elnino
# mydata = Octane.X
n = nrow(mydata)
p = ncol(mydata)

plot(as.numeric(mydata[1,]), col=grey(.6), type='l',
     ylim=(c(min(mydata), max(mydata)+.1)))
for(i in 2:n){
  lines(as.numeric(mydata[i,]), col=grey(.6))
}

# take a B-spline basis and orthonormalize it
B = bs(1:p, knots=seq(0,p-1,by=2), degree=3)
B1 = gramSchmidt(B)$Q

# transform data... diff between two bins is 1
Xt = as.matrix(mydata[,-1]) %*% B1[-1,]

# do depth PCA
source('misc_functions.R')
pcarank = PcaRank(Xt, k=2, proj=1e2)
Xthat = matrix(spatial.median(Xt, delta=1e-5)$mu, nrow=n, ncol=ncol(B1), byrow=T) + 
  pcarank@scores %*% t(pcarank@loadings)
Xhat = Xthat %*% t(B1)
Xhat[,1] = as.numeric(mydata[,1])

sdrank = pcarank@sd
odrank = pcarank@od
csdrank = pcarank@cutoff.sd
codrank = pcarank@cutoff.od
indices = 1:n
which.ind = which(sdrank > csdrank & odrank > codrank)

## distance-distance plots
pdf('Elnino_functional1.pdf',width=4, height=4)
plot(as.numeric(mydata[1,]), col=grey(.6), type='l',xaxt='n',
     ylab="Sea surface temperature", xlab="Month",
     ylim=(c(min(mydata), max(mydata)+.1)))
axis(side=1, at=1:12, labels=c('Jun','Jul','Aug','Sep','Oct','Nov','Dec',"Jan","Feb",'Mar','Apr','May'))
for(i in 2:n){
  lines(as.numeric(mydata[i,]), col=grey(.6))
}
for(i in which.ind){
  lines(as.numeric(mydata[i,]), col="black")
}
dev.off()

pdf('Elnino_functional2.pdf',width=4, height=4)
plot(as.numeric(Xhat[1,]), col=grey(.6), type='l',xaxt='n',
     ylab="Sea surface temperature", xlab="Month",
     ylim=(c(min(mydata), max(mydata)+.1)))
axis(side=1, at=1:12, labels=c('Jun','Jul','Aug','Sep','Oct','Nov','Dec',"Jan","Feb",'Mar','Apr','May'))
for(i in 2:n){
  lines(as.numeric(Xhat[i,]), col=grey(.6))
}
for(i in which.ind){
  lines(as.numeric(Xhat[i,]), col="black")
}
dev.off()

pdf('Elnino_functional3.pdf',width=4, height=4)
plot(sdrank, odrank,
     xlim=c(0,1.2*max(sdrank)), ylim=c(0,1.2*max(odrank)),
     pch=19, cex=.5,
     xlab="Score distance", ylab="Orthogonal distance")
abline(v=csdrank, col="red")
abline(h=codrank, col="red")
if(length(which.ind>0)){
  text(sdrank[which.ind]-.5, odrank[which.ind]-.5,
       paste0(indices[which.ind]+1949,"-",indices[which.ind]+50), pos=4, cex=.7)
}
dev.off()



