## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

set.seed(03092016)
df.list = c(5,6,10,15,25,Inf)
p.list = c(2,5,10,20)

ARE.mat1 = matrix(0, nrow=length(p.list), ncol=length(df.list))
for(i in 1:length(p.list)){
  for(j in 1:length(df.list)){
    pi = p.list[i]
    dfj = df.list[j]
    
    # ARE for PD
    X = matrix(rt(1e6*pi, df=dfj), ncol=pi)
    z = sqrt(rowSums(X^2))
    c = qnorm(.75)
    uz = z^2 / (c+z)^2
    upz = 2*c*z / (c+z)^3
    numerat = mean(pi*uz + upz*z)^2
    denom = mean(uz^2)
    ARE.mat1[i,j] = numerat/((pi)^2*denom)
  }
}

ARE.mat2 = matrix(0, nrow=length(p.list), ncol=length(df.list))
for(i in 1:length(p.list)){
  for(j in 1:length(df.list)){
    pi = p.list[i]
    dfj = df.list[j]
    
    # ARE for HSD
    X = matrix(rt(1e6*pi, df=dfj), ncol=pi)
    z = sqrt(rowSums(X^2))
    c = qnorm(.75)
    Fz = pt(z, df=dfj)
    uz = Fz^2
    upz = 2*Fz * dt(z, df=dfj)
    numerat = mean(z^4) * mean(pi*uz + upz*z)^2
    denom = mean(uz^2)
    ARE.mat2[i,j] = numerat/((pi^2+2*pi)^2*denom)
  }
}

cbind(t(ARE.mat1), t(ARE.mat2))

library(parallel)
library(doSNOW)

# Functions
ARE.norm = function(rho, depth, ns=1e4){
  require(fda.usc)
  Z = matrix(rnorm(2*ns), ncol=2)
  
  # get depth values
  if(depth=='HS'){
    DZ = mdepth.HS(Z,Z)$dep
  }
  else if(depth=='MhD'){
    DZ = mdepth.MhD(Z,Z)$dep
  }
  else {
    DZ = mdepth.RP(Z,Z)$dep
  }
  DZ = mdepth.RP(Z,Z)$dep
  DZ = max(DZ) - DZ
  
  
  DZ=1
  
  # calculate required quantities
  Z2 = Z*Z
  Z.Sig.Z = Z2[,1] + rho*Z2[,2]
  (E1 = mean(DZ^4 * Z2[,1] * Z2[,2] / Z.Sig.Z^2))
  (E2 = mean(DZ^2 * (Z2[,1] - rho*Z2[,2]) / Z.Sig.Z)^2)
  ARE = E2 / (E1*(1-rho)^2)
  
  ARE  
}

ARE.norm(.5, 'HS', ns=1e3)
ARE.norm(.25, 'MhD', ns=1e4)
ARE.norm(.5, 'RP', ns=1e4)

# Z = matrix(rnorm(2*1e4), ncol=2)
# r = sqrt(rowSums(Z*Z))
# mean(r^4)

lamS1 = mean(DZ^2*Z2[,1]/Z.Sig.Z)
lamS2 = rho*mean(DZ^2*Z2[,2]/Z.Sig.Z)
(AV1 = rho/(lamS1-lamS2)^2 * mean(DZ^4*Z2[,1]*Z2[,2]/Z.Sig.Z^2))
(AV2 = rho/(1-rho)^2)
AV2/AV1
