#Delta.list <- res[,r]$D
#Sigma.list <- res[,r]$S



stoppingCriterion0Para <- function(Tau.hat, Delta.list, Sigma.list, elements = NULL){
  
  d <- nrow(Tau.hat)
  l.ij.mat <- t(combn(1:d,2))
  p <- d * (d - 1) / 2
  
  tau.hat <- Tau.hat[l.ij.mat]
  
  
  ##################################
  cl <- makeCluster(detectCores()-1)
  t1 <- Sys.time()
  
  clusterExport(cl, c("Tau.hat", "Delta.list", "Sigma.list" , "tau.hat", "l.ij.mat", "d", "p"), envir = environment())
  clusterEvalQ(cl,library("Matrix"))
  # To use matrix.trace, should export the package matrixcalc also... or use matrixcalc::matrix.trace (?)
    
  t2 <- Sys.time()
  print(difftime(t2,t1))
  ##################################
  
  
  if(is.null(elements)){
    elements <- 1:d
  }
  
  t1 <- Sys.time()
  
  ones <- which(elements == 1)
  
  if(length(ones) > 0){
    elements.red <- elements[-ones]
  }else{
    elements.red <- elements
  }
  
  al <- parSapply(cl, elements.red, function(i){
    Delta <- Delta.list[[i]]
    Delta2 <- tcrossprod(Delta)
    
    clusters <- c(Delta %*% 1:dim(Delta)[2])
    ma <- max(clusters)
    
    double.index <- matrix(clusters[l.ij.mat], p, 2)
    double.index <- t(apply(X = double.index, MARGIN = 1, FUN = sort))
    
    l.blocks <- (double.index[,1] - 1) * ma - (double.index[,1] - 1) * (double.index[,1] - 2) / 2 + double.index[, 2] - double.index[, 1] + 1 
    
    B <- sparseMatrix(i = 1:p, j = l.blocks)
    if(ncol(B) > 1){
      B <- B[,unique(summary(B)$j)]
    }

    csB <- colSums(B)

    Gamma <- B %*% (t(B) / csB)
    
    #image(t(Gamma[p:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))
    #image(t((B %*% (t(B) / csB))[p:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))
    
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    vec <- crossprod(tau.tilde - tau.hat, sparseMatrix(1:p,1:p) - Gamma)
    
    R <- as.numeric(tcrossprod(vec,Sigma.list[[i]]^(-1) * vec))
    
    #print(matrix.trace(Gamma))
    #1 - pchisq(R, p - matrix.trace(Gamma))
    1 - pchisq(R, p - sum(diag(Gamma)))
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
