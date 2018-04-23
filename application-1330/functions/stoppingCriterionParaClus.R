#Delta.list <- res[,r]$D
#Sigma.list <- res[,r]$S



stoppingCriterionPara <- function(Tau.hat, clus.list, Sigma.list, elements = NULL, w = .5){
  
  d <- nrow(Tau.hat)
  l.ij.mat <- t(combn(1:d,2))
  p <- d * (d - 1) / 2
  
  tau.hat <- Tau.hat[l.ij.mat]
  
  if(is.null(elements)){
    elements <- 1:d
  }
  
  
  ##################################
  cl <- makeCluster(detectCores()-1)
  t1 <- Sys.time()
  
  ij.l.mat <- matrix(nrow=d,ncol=d)
  ij.l.mat[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(1:p,1:p)
  
  clusterExport(cl, c("Tau.hat", "clus.list", "Sigma.list" , "tau.hat", "l.ij.mat", "d", "p", "w"), envir = environment())
  clusterEvalQ(cl,library("Matrix"))
  # To use matrix.trace, should export the package matrixcalc also... or use matrixcalc::matrix.trace (?)
  
  t2 <- Sys.time()
  print(difftime(t2,t1))
  ##################################
  
  
  t1 <- Sys.time()
  
  ones <- which(elements == 1)
  
  if(length(ones) > 0){
    elements.red <- elements[-ones]
  }else{
    elements.red <- elements
  }
  
  al <- parSapply(cl, elements.red, function(i){
    clus <- clus.list[[i]]
    mem <- which(clus %in% pair)
    clus[mem] <- pair[1]
    clus <- frank(clus,ties.method="dense")
    ind <- do.call(what = rbind, lapply(seq_along(mem), function(i){
      cbind(rep(mem[i],d-i),(1:d)[-mem[1:i]])
    }))
    tau.tilde[ij.l.mat[ind]] <- Tt[cbind(clus[ind[,1]],clus[ind[,2]])]

    vec <- tau.tilde - tau.hat
    
    Si <- ((1-w)*Sigma.list[[i]] + w*mean(Sigma.list[[i]]))^(-1)
    R <- sum(vec^2*Si)
    
    #print(matrix.trace(Gamma))
    #1 - pchisq(R, p - matrix.trace(Gamma))
    K  <- d-i+1
    L <- K*(K-1)/2 + sum(sapply(clus, function(j){sum(clus==j) > 1}))
    1 - pchisq(R, p - L)
  })
  t2 <- Sys.time()
  print(difftime(t2,t1))
  
  stopCluster(cl)
  
  if(length(ones) > 0){
    new.al <- rep(1, length(elements))
    new.al[elements.red] <- al
    return(new.al)  
  }else{
    return(al)
  }
}

