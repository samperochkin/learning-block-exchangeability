stoppingCriterion <- function(Tau.hat, clus.list, Sigma.list, elements = NULL, w = .5){
  
  d <- nrow(Tau.hat)
  l.ij.mat <- t(combn(1:d,2))
  p <- d * (d - 1) / 2
  
  tau.hat <- Tau.hat[l.ij.mat]
  
  if(is.null(elements)){
    elements <- 1:d
  }
  
  
  ij.l.mat <- matrix(nrow=d,ncol=d)
  ij.l.mat[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(1:p,1:p)

  ones <- which(elements == 1)
  
  if(length(ones) > 0){
    elements.red <- elements[-ones]
  }else{
    elements.red <- elements
  }
  
  al <- sapply(elements.red, function(i){
    clus <- clus.list[[i]]
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
    vec <- tau.tilde - tau.hat
    
    Si <- ((1-w)*Sigma.list[[i]] + w*mean(Sigma.list[[i]]))^(-1)
    R <- sum(vec^2*Si)
    
    #print(matrix.trace(Gamma))
    #1 - pchisq(R, p - matrix.trace(Gamma))
    L <- K*(K-1)/2 + sum(sapply(cluss, length) > 1)
    1 - pchisq(R, p - L)
  })

  if(length(ones) > 0){
    new.al <- rep(1, length(elements))
    new.al[elements.red] <- al
    return(new.al)  
  }else{
    return(al)
  }
}

