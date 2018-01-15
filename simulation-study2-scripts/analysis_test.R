
num.rep <- 500
# load scenarios
source("Tau_generators.R")

n.scenario <- c(125,250,500)
w.vec <- c(0,.25,.5,.75,1)
scenario.grid <- expand.grid(list(n = n.scenario, Tau.num = 1:length(generator.list)))
scenario.grid2 <- expand.grid(list(w = w.vec, n = n.scenario, Tau.num = 1:length(generator.list)))

true.models <- c(17,17,13,13)

maxTrue.mat <- matrix(,0, num.rep)
bestSE.mat <- matrix(,0, num.rep)
hatSE.mat <- matrix(,0, num.rep)


for(i in 1:12){
  maxTrue.submat <- matrix(,length(w.vec), 0)
  bestSE.submat <- matrix(,length(w.vec), 0)
  hatSE.submat <- NULL
  
  Tau <- generator.list[[scenario.grid[i,2]]]
  
  for(j in 1:5){
    print(j)
    load(paste("res",i, "_", j, sep = ""))
    maxTrue.partial <- sapply(1:100, function(jj){
      sapply(1:length(w.vec), function(k){
        which(sapply(res[[jj]][,k]$Delta, function(x){sum(true.delta[[scenario.grid[i,2]]] - tcrossprod(x) < -.5) > 0}))[1] - 1
      })
    })
    bestSE.partial <- sapply(1:100, function(jj){
      sapply(1:length(w.vec), function(k){
        min(sapply(res[[jj]][,k]$Tau, function(x){sum((x - Tau)^2)}))
      })
    })
    hatSE.partial <- sapply(1:100, function(jj){
      sum((res[[jj]][,1]$Tau[[1]] - Tau)^2)
    })
    
    maxTrue.submat <- cbind(maxTrue.submat,maxTrue.partial)
    bestSE.submat <- cbind(bestSE.submat,bestSE.partial)
    hatSE.submat <- c(hatSE.submat, hatSE.partial)
  }
  
  maxTrue.mat <- rbind(maxTrue.mat,maxTrue.submat)
  bestSE.mat <- rbind(bestSE.mat,bestSE.submat)
  hatSE.mat <- rbind(hatSE.mat, hatSE.submat)
}

vec <- c(rep(17,nrow(scenario.grid2)/2),rep(13,nrow(scenario.grid2)/2))
getTrue.mat <- maxTrue.mat == vec
save(list = c("getTrue.mat","maxTrue.mat","bestSE.mat","hatSE.mat"), file = "perf1_mats.R")
