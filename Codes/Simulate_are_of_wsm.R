## calculate asymptotic relative efficiency of weighted spatial median
setwd("C:/Study/My projects/Depth-scatter/Codes")
rm(list=ls());
source('misc_functions.R')

ARE.wsm = function(p, dist, df=NULL, rho=0){
  # covariance matrix
  Sigma = diag(p)
  for(i in 2:p){
    for(j in 1:(i-1)){
      Sigma[i,j] = rho^abs(i-j)
    }
  }
  Sigma = Sigma + t(Sigma) - diag(p)
  
  # generate random entries
  if(dist=="norm"){
    e.mat = my.mvrnorm(1e4, mu=rep(0,p), Sigma=Sigma)
  } else{
    e.mat = my.mvrt(1e4, df=df, mu=rep(0,p), Sigma=Sigma)
  }
  A.list = list()
  B.list = list()
  Aw.list = list()
  Bw.list = list()
  
  for(i in 1:1e4){
    z = e.mat[i,]
    normz = sqrt(sum(z^2))
    ztz = z %*% t(z)
    wz = 1/(normz)
    A.list[[i]] = (diag(p) - ztz/normz^2) / normz
    B.list[[i]] = ztz / normz^2
    Aw.list[[i]] = (diag(p) - ztz/normz^2) * wz / normz
    Bw.list[[i]] = wz^2 * ztz / normz^2
  }
  
  # calculate quantities
  A = apply(simplify2array(A.list), 1:2, mean)
  B = apply(simplify2array(B.list), 1:2, mean)
  Aw = apply(simplify2array(Aw.list), 1:2, mean)
  Bw = apply(simplify2array(Bw.list), 1:2, mean)
  
  (det(solve(A) %*% B %*% solve(A)) / det(solve(Aw) %*% Bw %*% solve(Aw)))^(1/p)
}

set.seed(8202016)
c(ARE.wsm(5,"norm", 0),ARE.wsm(5,"t",df=3, 0),
  ARE.wsm(5,"t", df=5, 0),ARE.wsm(5,"t", df=10, 0),ARE.wsm(5,"t", df=20, 0))
c(ARE.wsm(10,"norm", 0),ARE.wsm(10,"t",df=3, 0),
  ARE.wsm(10,"t", df=5, 0),ARE.wsm(10,"t", df=10, 0),ARE.wsm(10,"t", df=20, 0))
c(ARE.wsm(20,"norm", 0),ARE.wsm(20,"t",df=3, 0),
  ARE.wsm(20,"t", df=5, 0),ARE.wsm(20,"t", df=10, 0),ARE.wsm(20,"t", df=20, 0))
c(ARE.wsm(50,"norm", 0),ARE.wsm(50,"t",df=3, 0),
  ARE.wsm(50,"t", df=5, 0),ARE.wsm(50,"t", df=10, 0),ARE.wsm(50,"t", df=20, 0))
