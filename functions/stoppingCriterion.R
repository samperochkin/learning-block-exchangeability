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
