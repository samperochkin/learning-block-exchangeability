computeTt <- function(new.Delta, Theta.hat, l.ij.mat){
  
  d <- dim(new.Delta)[1]
  p <- d * (d-1) / 2
  K <- dim(new.Delta)[2]
  
  new.Theta <- Theta.hat
  
  clusters <- c(tcrossprod(1:K,new.Delta))
  clusters <- lapply(1:K, function(k){
    which(clusters == k)
  })
  
  for(k in 1:length(clusters)){
    new.Theta[clusters[[k]]] <- mean(Theta.hat[clusters[[k]]])
  }
  
  new.Theta
}