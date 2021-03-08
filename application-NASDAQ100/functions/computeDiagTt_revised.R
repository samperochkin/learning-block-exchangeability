computeTt <- function(new.Delta, theta.hat, l.ij.mat){
  
  d <- dim(new.Delta)[1]
  p <- d * (d-1) / 2
  K <- dim(new.Delta)[2]
  
  new.theta <- theta.hat
  
  
  clusters <- c(tcrossprod(1:K,new.Delta))
  clus2 <- matrix(clusters[c(l.ij.mat)], ncol = 2)
  clus2 <- t(apply(clus2, 1, sort))
  clus1 <- clus2[,1]*(K+1) + clus2[,2]
  uclus1 <- unique(clus1)
  
  for(k in uclus1){
    new.theta[clus1 == k] <- mean(theta.hat[clus1 == k])
  }
  
  new.theta
}

#image(t(new.Theta[p:1,]), zlim = c(min(new.Theta[p:1,]),max(new.Theta[p:1,])), axes=FALSE, col=colfunc(100))
#identical(new.Theta, Theta.list[[3]])
