## Sim_efficiency: Asymptotic relative efficiency for depth-based scatter estimator
set.seed(03092016)
df.list = c(5,6,10,15,25,Inf)
p.list = c(2,5,10,20)
N = 1e6

ARE.mat1 = matrix(0, nrow=length(p.list), ncol=length(df.list))
for(i in 1:length(p.list)){
    pi = p.list[i]
    lambda = rev(1:pi)
    sqrtLambda = sqrt(lambda)
    # lamratio = outer(lambda,lambda,"*")/outer(lambda,lambda,"-")^2
    # diag(lamratio) = 0
    for(j in 1:length(df.list)){
        dfj = df.list[j]
        
        # calculate MAD
        if(dfj==Inf){
            c = qnorm(.75)
        } else{
            c = 2*sqrt(dfj)/((dfj-1)*beta(dfj/2,1/2))
        }
        
        # ARE for PD
        Z = matrix(rt(N*pi, df=dfj), ncol=pi)
        U = sweep(Z, 1, nZ, FUN="/")
        nZ = sqrt(rowSums(Z^2))
        W = nZ/(1+nZ/c)
        X = Z %*% diag(sqrtLambda)
        S = sweep(X, 1, sqrt(rowSums(X^2)), FUN="/")
        WS = sweep(S, 1, W, FUN="*")
        lamw = diag(cov(WS))
        lamwdiff = outer(lamw, lamw, FUN = "-")
        W2SS1 = sweep(WS, 1, WS[,1], FUN="*")
        d.vec = colMeans(W2SS1^2)[-1]/lamwdiff[-1,1]^2
        lam = diag(cov(S))
        lamdiff = outer(lam, lam, FUN = "-")
        SS1 = sweep(S, 1, S[,1], FUN="*")
        n.vec = colMeans(SS1^2)[-1]/lamdiff[-1,1]^2
        ARE.mat1[i,j] = sum(n.vec) / sum(d.vec)
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
