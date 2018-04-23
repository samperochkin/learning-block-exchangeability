criterion <- function(clus, w = .5, Tau.hat, Sigma.hat, l.ij.mat, ij.l.mat){
  
  K <- max(clus)
  cluss <- lapply(1:K, function(k){
    which(clus==k)
  })

  Tt <- matrix(sapply(1:K, function(k){
    sapply(1:K, function(kk){
      if(k == kk){
        ll <- length(cluss[[k]])
        sum(Tau.hat[cluss[[k]],cluss[[k]]] - diag(ll))/(ll*(ll-1))
      }else{
        mean(Tau.hat[cluss[[k]],cluss[[kk]]])
      }
    })
  }),K,K)
    

  tau.tilde <- Tt[cbind(clus[l.ij.mat[,1]],clus[l.ij.mat[,2]])]
  vec <- tau.tilde - Tau.hat[l.ij.mat]
  
  # average Sigma and apply criterion!
  
  R <- sum(vec^2/((1-w)*Sigma.hat + w*mean(Sigma.hat)))
  L <- K*(K-1)/2 + sum(sapply(cluss, length) > 1)
  1 - pchisq(R, nrow(l.ij.mat) - L)
}

