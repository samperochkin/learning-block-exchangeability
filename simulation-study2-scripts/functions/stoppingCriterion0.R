stoppingCriterion0 <- function(Tau.hat, Delta.list, Sigma.list){
  
  d <- nrow(Tau.hat)
  l.ij.mat <- t(combn(1:d,2))
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
    1 - pchisq(R, p - matrix.trace(Gamma))
  })
  
  c(1,al)
}