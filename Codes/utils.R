## utils: utility functions for ease of use

# function giving matrix of ones
ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

## computes the spatial median
spatial.median <- function(x, delta=1e-3)
{
  dime = dim(x)
  n=dime[1]
  p=dime[2]
  delta1=delta*sqrt(p)
  mu0=apply(x,2,median)
  h=delta1+1
  tt=0
  while(h>delta1)
  {
    tt=tt+1
    TT=matrix(mu0,n,p,byrow=TRUE)
    U=(x-TT)^2
    w=sqrt(apply(U,1,sum))
    w0=median(w)
    ep=delta*w0
    
    z=(w<=ep)
    w[z]=ep
    w[!z]=1/w[!z]
    w=w/sum(w)
    x1=x
    for(i in 1:n)
      x1[i,]=w[i]*x[i,]
    mu=apply(x1,2,sum)
    h=sqrt(sum((mu-mu0)^2))
    mu0=mu
  }
  out=list(mu=mu0,ep=ep)
  out
}

# compute weighted spatial medians
depth.median <- function(x, d, delta=1e-3)
{
  dime = dim(x)
  n=dime[1]
  p=dime[2]
  delta1=delta*sqrt(p)
  mu0=apply(x,2,median)
  h=delta1+1
  tt=0
  
  while(h>delta1)
  {
    tt=tt+1
    TT=matrix(mu0,n,p,byrow=TRUE)
    U=(x-TT)^2
    w=sqrt(apply(U,1,sum))
    w0=median(w)
    ep=delta*w0
    
    z=(w<=ep)
    w[z]=ep
    w[!z]=1/w[!z] * d
    w=w/sum(w)
    x1=x
    for(i in 1:n)
      x1[i,]=w[i]*x[i,]
    mu=apply(x1,2,sum)
    h=sqrt(sum((mu-mu0)^2))
    mu0=mu
  }
  out=list(mu=mu0,ep=ep)
  out
}


# compute Tyler's shape matrix, optionally with weights
TylerSig = function(X, tol=1e-5, maxit=100, weight=NULL){
    n = nrow(X); p = ncol(X)
    iSig = diag(rep(1,p))
    
    # whether to use depth weights
    if(is.null(weight)){
        weight = rep(1,n)
    }
    
    for(i in 1:maxit){
        iiSig = matrix(0,p,p)
        inv.iSig = solve(iSig)
        for(j in 1:n){
            xj = as.matrix(X[j,])
            iiSig = iiSig + weight[j]^2 * (xj %*% t(xj))/ as.numeric(t(xj) %*% inv.iSig %*% xj)
        }
        
        iiSig = iiSig/det(iiSig)^(1/p)
        if(norm(iSig - iiSig, type="F") < tol){
            break
        }
        else{
            iSig = iiSig
        }
    }
    
    iSig
}