## Sim_efficiency: Finite sample efficiency for depth-based scatter estimator
## comparison with sign covariance matrix (SCM)

setwd("C:/Users/smajumdar/Documents/Depth-scatter1/Codes")
Required.Packages <- c("mvtnorm","parallel","doSNOW", "fastM","ddalpha","data.table")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})
source('misc_functions.R')

n = c(20,50,100,300,500)
p = 4
iter = 1e3
i=3
ncores=6

## Functions
add_contamination = function(X, eps=eps, M=1e3){
    N = nrow(X)
    Nc = ceiling(N*eps)
    Xnew = X
    if(Nc>0){
        set.seed(Nc)
        eps_rows = sample(1:N, Nc)
        for(row in eps_rows){
            Xnew[row,1] = X[row,1] + max(abs(X))*M
        }
    }
    Xnew
}

## simulate for normal dist
FSE.sim = function(n, p, eps=0, iter=1e3, dist="norm", df=NULL, ncores=detectCores()){
    set.seed(12182014)
    v = c(rep(0,p-1), 1)
    lam = 1:p
    Sigma = diag(lam)
    nn = length(n)
    ones = rep(1,p)
    MSE.mat = matrix(0, nrow=nn, ncol=9)
    time.mat = matrix(0, nrow=nn, ncol=9)
    
    # generate data for largest n
    X_list = lapply(1:iter, function(j){
        set.seed(1e3*j)
        if(dist=="norm"){
            matrix(rnorm(p*max(n)), ncol=p) %*% sqrt(Sigma)
        } else{
            rmvt(max(n), sigma=Sigma, df=df)
        }
    })
    
    for(i in 1:nn){
        
        ni = n[i]
        cat("Doing n =",ni,"-> ")
        # create data matrices, add contamination if needed
        iX_list = lapply(X_list, function(x) add_contamination(x[1:ni,], eps=eps))
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
        MSE.mat[i,] = MSE.vec
        time.mat[i,] = it
        cat("done\n")
    }
    
    list(MSE.mat, time.mat)
}

# For p=4
n.vec = seq(50,500,by=50)
norm.table <- FSE.sim(n=n.vec, p=4, iter=1e3, ncores=7)
t3.table <- FSE.sim(n=n.vec, p=4, dist="t", df=3, iter=1e3, ncores=7)
t10.table <- FSE.sim(n=n.vec, p=4, dist="t", df=10, iter=1e3, ncores=7)
t20.table <- FSE.sim(n=n.vec, p=4, dist="t", df=20, iter=1e3, ncores=7)
saveRDS(list(norm.table, t3.table, t10.table, t20.table), file="normt4.rds")

# contamination added
c1.table = FSE.sim(n=n.vec, p=4, eps=.1, iter=1e3, ncores=7)
c3.table = FSE.sim(n=n.vec, p=4, eps=.3, iter=1e3, ncores=7)
saveRDS(list(c1.table, c3.table), file="normcont4.rds")
