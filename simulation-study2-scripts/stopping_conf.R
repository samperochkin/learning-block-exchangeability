# stopping loop
stoppingCriterion <- function(Tau.hat, Delta.list, Sigma.list, alpha, l.ij.mat){
  
  d <- dim(Tau.hat)[1]
  p <- d * (d - 1) / 2
  
  tau.hat <- Tau.hat[l.ij.mat]
  
  al <- sapply(2:d, function(i){
    #print(i)
    Delta <- Delta.list[[i]]
    Delta2 <- tcrossprod(Delta)
    
    clusters <- c(Delta %*% 1:dim(Delta)[2])
    ma <- max(clusters)
    
    double.index <- matrix(clusters[l.ij.mat], p, 2)
    double.index <- t(apply(X = double.index, MARGIN = 1, FUN = sort))
    
    l.blocks <- (double.index[,1] - 1) * ma - (double.index[,1] - 1) * (double.index[,1] - 2) / 2 + double.index[, 2] - double.index[, 1] + 1 
    
    B <- matrix(0,p,max(l.blocks))
    B[cbind(1:p,l.blocks)] <- 1
    
    Gamma <- B%*%t(B)%*%ginv(B%*%t(B))
    
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    vec <- tau.tilde - tau.hat
    R <- crossprod(vec,diag(p) - Gamma)  %*% solve(Sigma.list[[i]]) %*% vec
    
    #print(matrix.trace(Gamma))
    1 - pchisq(R, p - matrix.trace(Gamma))
  })
  
  #sapply(alpha, function(a){rev(which(cumsum(c(1,al > a))/(1:d) == 1))[1]})
  #sapply(alpha, function(a){rev(which(c(1,al) > a))[1]})
  c(1,al)
}


num.rep <- 500
# load scenarios
#source("C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/R code/Simulation_study/Full_study/Full_study_project/Tau_generators.R")
source("Tau_generators.R")



n.scenario <- c(125,250,500)
w.vec <- c(0,.25,.5,.75,1)

scenario.grid <- expand.grid(list(n = n.scenario, Tau.num = 1:length(generator.list)))
scenario.grid2 <- expand.grid(list(w = w.vec, n = n.scenario, Tau.num = 1:length(generator.list)))

true.models <- c(17,17,13,13)

l.ij.fun <- function(l,d){
  i <- d-sum(l<=cumsum(d-(1:(d-1))))
  j <- l - ((i-1)*d - i*(i-1)/2 - i)
  c(i,j)
}

library(parallel)
cl <- makeCluster(detectCores()-1)

# export libraries
clusterEvalQ(cl,library("MASS"))
clusterEvalQ(cl,library("matrixcalc"))
clusterEvalQ(cl,library("pcaPP"))
clusterEvalQ(cl,library("data.table"))

# export objects
clusterExport(cl, varlist = c("stoppingCriterion", "l.ij.fun"))

parSapply(cl, 1:12, function(i){
#for(i in 1:12){
  alpha1.mat <- matrix(, 20,0)
  alpha2.mat <- matrix(, 20,0)
  alpha3.mat <- matrix(, 20,0)

  for(j in 1:5){
    #print(j)
    load(paste("res",i, "_", j, sep = ""))
    
    d <- length(res[[1]][,1]$Delta)
    p <- d * (d-1) / 2
    l.ij.mat <- t(sapply(1:p,function(l){return(l.ij.fun(l, d))}))
    
    alpha1.partial <- matrix(0,20,100)
    alpha2.partial <- matrix(0,20,100)
    
    alpha1.partial <- sapply(1:100, function(jj){
      #print(jj)
      stoppingCriterion(res[[jj]][,4]$Tau[[1]], res[[jj]][,4]$Delta, res[[jj]][,4]$Sigma, 0, l.ij.mat)
    })

    alpha2.partial <- sapply(1:100, function(jj){
      #print(jj)
      SS <- res[[jj]][,4]$Sigma
      SS <- sapply(SS, function(S){.75 * diag(diag(S)) + .25 * S}, simplify = FALSE)
      stoppingCriterion(res[[jj]][,4]$Tau[[1]], res[[jj]][,4]$Delta, SS, 0, l.ij.mat)
    })
    
    alpha1.mat <- cbind(alpha1.mat, alpha1.partial)
    alpha2.mat <- cbind(alpha2.mat, alpha2.partial)
  }
  
  saveRDS(alpha1.mat, paste0("alpha0_",i))
  saveRDS(alpha2.mat, paste0("alpha75_",i))
})

stopCluster(cl)
