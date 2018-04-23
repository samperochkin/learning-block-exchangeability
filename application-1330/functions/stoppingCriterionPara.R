#Delta.list <- res[,r]$D
#Sigma.list <- res[,r]$S



stoppingCriterionPara <- function(Tau.hat, Delta.list, Sigma.list, elements = NULL, w = .5){
  
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
  
  clusterExport(cl, c("Tau.hat", "Delta.list", "Sigma.list" , "tau.hat", "l.ij.mat", "d", "p", "w"), envir = environment())
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
    Delta <- Delta.list[[i]]
    Delta2 <- tcrossprod(Delta)
    
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    vec <- tau.tilde - tau.hat
    
    Si <- ((1-w)*Sigma.list[[i]] + w*mean(Sigma.list[[i]]))^(-1)
    R <- as.numeric(crossprod(vec^2, Si))

    #print(matrix.trace(Gamma))
    #1 - pchisq(R, p - matrix.trace(Gamma))
    K  <- ncol(Delta)
    L <- K*(K-1)/2 + sum(colSums(Delta) > 1)
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
