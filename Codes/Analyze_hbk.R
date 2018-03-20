## Analyze_octane: analysis of octane data
rm(list=ls())
setwd("C:/Study/My projects/Depth-scatter/Codes")
source('misc_functions.R')

## Octane data
library(ellipse)
library(rrcov)
library(fda.usc)
set.seed(04112015)

## Analyze bus data
data(bus)
bus1 <- bus[, -9]
madbus <- apply(bus1, 2, mad)
bus2 <- sweep(bus1, 2, madbus, "/", check.margin = FALSE)
# bus2 = scale(bus2)

system.time(pca <- PcaClassic(bus2))
system.time(rpca <- PcaLocantore(bus2))
system.time(pcaHubert <- PcaHubert(bus2, k=17, kmax=17, mcd=FALSE))
system.time(pcamcd <- PcaCov(bus2, cov.control=CovControlMcd()))
#pcaogk <- PcaCov(bus2, cov.control=CovControlOgk())
system.time(pcarank <- PcaRank(bus2))

ev <- getEigenvalues(pca)
evrob <- getEigenvalues(rpca)
evhub <- getEigenvalues(pcaHubert)
evmcd <- getEigenvalues(pcamcd)
#evogk <- getEigenvalues(pcaogk)
evrank = getEigenvalues(pcarank)

uvar <- matrix(nrow=6, ncol=6)
svar <- sum(ev)
svarrob <- sum(evrob)
svarhub <- sum(evhub)
svarmcd <- sum(evmcd)
#svarogk <- sum(evogk)
svarrank = sum(evrank)

for(i in 1:6){
  uvar[i,1] <- i
  uvar[i,2] <- round((svar - sum(ev[1:i]))/svar, 3)
  uvar[i,3] <- round((svarrob - sum(evrob[1:i]))/svarrob, 3)
  uvar[i,4] <- round((svarhub - sum(evhub[1:i]))/svarhub, 3)
  uvar[i,5] <- round((svarmcd - sum(evmcd[1:i]))/svarmcd, 3)
#  uvar[i,6] <- round((svarogk - sum(evogk[1:i]))/svarogk, 3)
  uvar[i,6] <- round((svarrank - sum(evrank[1:i]))/svarrank, 3)
}
uvar <- as.data.frame(uvar)
names(uvar) <- c("q", "Classical","Spherical", "Hubert", "MCD", "Depth")
cat("\nBus data: proportion of unexplained variability for q components\n")
print(uvar)

## Reproduce Table 6.4 from Maronna et al. (2006), page 214, adding DCM
##
## Compute classical and robust PCA extracting only the first 3 components
## and take the squared orthogonal distances to the 3-dimensional hyperplane
##
pca3 <- PcaClassic(bus2, k=3)               # classical
rpca3 <- PcaLocantore(bus2, k=3)            # spherical (Locantore, 1999)
hpca3 <- PcaHubert(bus2, k=3)               # Hubert
mpca3 = PcaCov(bus2, cov.control=CovControlMcd(), k=3) # MCD
dpca3 = PcaRank(bus2, k=3)                  # Depth

dist <- pca3@od^2
rdist <- rpca3@od^2
hdist <- hpca3@od^2
mdist = mpca3@od^2
ddist = dpca3@od^2

## calculate the quantiles of the distances to the 3-dimensional hyperplane
qclass  <- round(quantile(dist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
qspc <- round(quantile(rdist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
qhubert <- round(quantile(hdist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
qmcd <- round(quantile(mdist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)
qdepth <- round(quantile(ddist, probs = seq(0, 1, 0.1)[-c(1,11)]), 1)

qq <- cbind(rbind(qclass, qspc, qhubert, qmcd, qdepth),
            round(c(max(dist), max(rdist), max(hdist), max(mdist), max(ddist)), 0))
colnames(qq)[10] <- "Max"
rownames(qq) <- c("Classical", "Spherical", "Hubert", "MCD", "Depth")
cat("\nBus data: quantiles of distances to hyperplane\n")
print(qq)
