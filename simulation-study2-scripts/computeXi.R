final.mat <- do.call(rbind,sapply(1:12, function(i){
  alphas <- readRDS(paste0("alpha0_",i))
  
  sapply(1:ncol(alphas), function(k){
    sapply(c(.05,.1,.25,.5), function(aa){
      rev(which(alphas[,k] > aa))[1]
    })
  })
}, simplify = FALSE))

library(parallel)
cl <- makeCluster(detectCores()-1)
clusterExport(cl, "final.mat")

xi.mat <- t(parSapply(cl,1:48, function(ii){
  i <- ceiling(ii/4)
  source("Tau_generators.R")
  Tau <- generator.list[[ceiling(ii/12)]]
  xi.sub <- NULL
  for(j in 1:5){
    load(paste0("res",i,"_",j))
    xi.partial <- sapply(1:100, function(jj){
      ff <- final.mat[ii,(j-1) * 100 + jj]
      1 - sum((res[[jj]][,4]$Tau[[ff]] - Tau)^2)/sum((res[[jj]][,4]$Tau[[1]] - Tau)^2)
    })
    xi.sub <- c(xi.sub,xi.partial)
  }
  xi.sub
}))

stopCluster(cl)

save(list = "xi.mat", file = "xi_mat.R")

