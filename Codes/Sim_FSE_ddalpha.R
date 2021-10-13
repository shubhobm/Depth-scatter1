## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

setwd("C:/Users/smajumdar/Documents/Depth-scatter1/Codes")
Required.Packages <- c("parallel","doSNOW", "fastM","ddalpha","data.table")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source('misc_functions.R')

n = c(20,50,100,300,500)
p = 4
iter = 1e3
i=3
ncores=6

## Functions

## simulate for normal dist
FSE.norm = function(n, p, iter=1e3, ncores=detectCores()){
    set.seed(12182014)
    v = c(rep(0,p-1), 1)
    lam = 1:p
    Sigma = diag(lam)
    nn = length(n)
    ones = rep(1,p)
    MSE.mat = matrix(0, nrow=nn, ncol=8)
    time.mat = matrix(0, nrow=nn, ncol=8)
    
    # generate data for largest n
    X_list = lapply(1:iter, function(j){
        set.seed(1e3*j)
        matrix(rnorm(p*max(n)), ncol=p) %*% sqrt(Sigma)
    })
    
    for(i in 1:nn){
        
        ni = n[i]
        cat("Doing n =",ni,"-> ")
        iX_list = lapply(X_list, function(x) x[1:ni,])
        iv_list = vector("list",9)
        it = rep(0,9)
        
        ## calculate all FSE values
        # vanilla PCA
        it[1] = system.time(iP_list <- lapply(iX_list, princomp))[3]
        iv_list[[1]] = sapply(iP_list, function(x) abs(sum(v * x$loadings[,1])))
        
        # SCM
        cl = makeCluster(ncores)
        registerDoSNOW(cl)
        scm_fun = function(j){
            iX = iX_list[[j]]
            iXnorm = sqrt(iX^2 %*% ones)
            iS = iX / (iXnorm %*% ones)
            iPsign = princomp(iS)
            list(iS, iPsign)
        }
        it[2] = system.time(scm_list <- foreach(j=1:iter) %dopar% scm_fun(j))[3]
        iS_list = lapply(scm_list, function(x) x[[1]])
        iv_list[[2]] = sapply(scm_list, function(x) abs(sum(v * x[[2]]$loadings[,1])))
        
        # Tyler's scatter
        tyl_fun = function(j){
            require(fastM)
            eigen(TYLERshape(iX_list[[j]])$Sigma)
        }
        it[3] = system.time(tyl_list <- foreach(j=1:iter) %dopar% tyl_fun(j))[3]
        iv_list[[3]] = sapply(tyl_list, function(x) abs(sum(v * x$vectors[,1])))
        
        # Halfspace depth
        dep_fun = function(j, depth=c("h","m","p")){
            require(ddalpha)
            iX = iX_list[[j]]
            if(depth=="h"){
                idep = depth.halfspace(iX,iX)
            } else if(depth=="m"){
                idep = depth.Mahalanobis(iX,iX)
            } else{
                idep = depth.projection(iX,iX)
            }
            idep = max(idep) - idep
            iXd = iS_list[[j]] * idep
            list(idep, princomp(iXd))
        }
        it[4] = system.time({
            dep_list <- foreach(j=1:iter) %dopar% dep_fun(j, depth="h")
        })[3]
        iv_list[[4]] = sapply(dep_list, function(x) abs(sum(v * x[[2]]$loadings[,1])))
        
        idepTyler = function(j){
            source('misc_functions.R')
            eigen(TylerSig(iX_list[[j]], weight=dep_list[[j]][[1]]))
        }
        it[5] = system.time(dep_list2 <- foreach(j=1:iter) %dopar% idepTyler(j))[3]
        iv_list[[5]] = sapply(dep_list2, function(x) abs(sum(v * x$vectors[,1])))
        
        # mahalanobis depth
        it[6] = system.time({
            dep_list <- foreach(j=1:iter) %dopar% dep_fun(j, depth="m")
        })[3]
        iv_list[[6]] = sapply(dep_list, function(x) abs(sum(v * x[[2]]$loadings[,1])))
        
        it[7] = system.time(dep_list2 <- foreach(j=1:iter) %dopar% idepTyler(j))[3]
        iv_list[[7]] = sapply(dep_list2, function(x) abs(sum(v * x$vectors[,1])))
        
        # projection depth
        it[8] = system.time({
            dep_list <- foreach(j=1:iter) %dopar% dep_fun(j, depth="p")
        })[3]
        iv_list[[8]] = sapply(dep_list, function(x) abs(sum(v * x[[2]]$loadings[,1])))
        
        it[9] = system.time(dep_list2 <- foreach(j=1:iter) %dopar% idepTyler(j))[3]
        iv_list[[9]] = sapply(dep_list2, function(x) abs(sum(v * x$vectors[,1])))
        
        stopCluster(cl)
        
        # get MSE and return
        eff.v = matrix(unlist(iv_list), ncol=9, byrow=F)
        (MSE.vec = apply(eff.v, 2, function(x) mean(acos(x)^2)))
        MSE.mat[i,] = MSE.vec[1]/MSE.vec[-1]
        time.mat[i,] = it[-1]
        cat("done\n")
    }
    
    list(MSE.mat, time.mat)
}

## simulate for t-distn
FSE.t = function(n, p, df, iter=1e3, ncores=detectCores()){
  set.seed(12182014)
  v = c(rep(0,p-1), 1)
  lam = 1:p
  Sigma = diag(lam)
  
  MSE.mat = matrix(0, nrow=length(n), ncol=8)
  for(i in 1:length(n)){
    
    # function to compute stuff 1000 times for a given n
    loopfun = function(j){
      source('misc_functions.R')
      require(fastM)
      require(fda.usc)
      require(mvtnorm)
      iv = rep(0,9)
      
      # get sample and construct sign matrix
      iX = rmvt(n[i], sigma=Sigma, df=df)
      iXnorm = sqrt(iX^2 %*% rep(1,p))
      iS = iX / (iXnorm %*% rep(1,p))
      
      # PCA on original sample
      iP = princomp(iX)
      iv[1] = abs(sum(v * iP$loadings[,1]))
      
      # PCA on SCM
      iPsign = princomp(iS)
      iv[2] = abs(sum(v * iPsign$loadings[,1]))
      
      # PCA on Tyler's cov matrix
      T = TYLERshape(iX)$Sigma
      iv[3] = abs(sum(v * eigen(T)$vectors[,1]))
      
      # PCA on DCM and depth-weighted tyler's scatter
      # Tukey's depth
      idep = mdepth.HS(iX,iX)$dep
      #idep = EPQD(iX,iX)[,p+1]
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[4] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[5] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      # Mahalanobis depth
      idep = mdepth.MhD(iX,iX)$dep
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[6] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[7] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      # Projection depth
      idep = mdepth.RP(iX,iX)$dep
      idep = max(idep) - idep
      iXd = iS * idep
      iPdepth = princomp(iXd)
      iv[8] = abs(sum(v * iPdepth$loadings[,1]))
      
      Td = TylerSig(iX, weight=idep)
      iv[9] = abs(sum(v * eigen(Td)$vectors[,1]))
      
      iv
    }
    
    # parallel code: compute MSE elements iter times
    cl = makeCluster(ncores)
    registerDoSNOW(cl)
    system.time(eff.v <- foreach(j=1:iter) %dopar% loopfun(j))
    stopCluster(cl)
    
    # get MSE and return
    eff.v = matrix(unlist(eff.v), ncol=9, byrow=T)
    (MSE.vec = apply(eff.v, 2, function(x) mean(acos(x)^2)))
    MSE.mat[i,] = MSE.vec[1]/MSE.vec[-1]
  }
  
  MSE.mat
}

# For p=4
n.vec = seq(50,500,by=50)
norm.table <- FSE.norm(n.vec, 4, 1e3, ncores=8)
saveRDS(norm.table, file="norm4.rds")
