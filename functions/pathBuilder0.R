pathBuilder0 <- function(Tau.hat, Theta.hat, n){
  
  #### fixed parameters ####
  
  d <- dim(Tau.hat)[1]
  p <- d * (d - 1) / 2
  
  #### vectorizing ####
  
  l.ij.mat <- t(sapply(1:p,function(l){
    i <- d-sum(l<=cumsum(d-(1:(d-1))))
    j <- l - ((i-1)*d - i*(i-1)/2 - i)
    c(i,j)
  }))
  
  m.lr.mat <- matrix(rep(1:p,2),p,2)
  
  # lazy way to avoid errors in previously written code
  lr.m.mat <- diag(1:p)

  m.shared.mat <- t(sapply(1:dim(m.lr.mat)[1], function(m){
    shared <- intersect(l.ij.mat[m.lr.mat[m, 1],],l.ij.mat[m.lr.mat[m, 2],])
    if(length(shared) == 0){
      return(c(0,0))
    }else if(length(shared) == 1){
      return(c(0,shared))
    }else{
      shared <- sort(shared)
      return(shared)
    }
  }))
  
  #### Additional initial values ####
  
  Tau.list <- list()
  Sigma.list <- list()
  Delta.list <- list()
  
  Tau.list[[1]] <- Tau.hat
  tau.hat <- Tau.hat[l.ij.mat]
  
  Sigma.hat <- buildSigma(Theta.hat, tau.hat, n)
  Sigma.list[[1]] <- Sigma.hat

  S <- Sigma.list[[1]]
  Si <- S^(-1)
  
  Delta <- diag(d)
  Delta.list[[1]] <- Delta
  
  
  #### Functions for the iterations ####
  
  DeltaCandidates <- function(Delta){
    K <- dim(Delta)[2]
    L <- (K*(K-1)/2)
    ij.pairs <- t(sapply(1:L,function(l){
      i <- K-sum(l<=cumsum(K-(1:(K-1))))
      j <- l - ((i-1)*K - i*(i-1)/2 - i)
      c(i,j)
    }))
    
    sapply(1:L, function(l){
      pair <- ij.pairs[l,]
      Delta[,pair[1]] <- rowSums(Delta[,pair])
      Delta[,-pair[2]]
    }, simplify = FALSE)
  }
  
  loss <- function(Delta, Tau.hat, Si, l.ij.mat, d){
    # Avoiding using Moore-Penrose for real.
    Delta2 <- tcrossprod(Delta)
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    vec <- tau.hat - tau.tilde
    crossprod(vec^2, Si)
  }
  
  #### Iterations ####
  # i <- 3
  for(i in 2:(d-1)){
    Delta.candidates <- DeltaCandidates(Delta.list[[i-1]])
    
    new.Delta <- Delta.candidates[[which.min(sapply(Delta.candidates, function(Delta){loss(Delta, Tau.hat, Si, l.ij.mat, d)}))]]
    
    t1 <- Sys.time()
    Delta2 <- tcrossprod(new.Delta)
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    Theta.tilde <- computeTt(new.Delta, Theta.hat, m.lr.mat, lr.m.mat, m.shared.mat, l.ij.mat)
    Sigma.list[[i]] <- buildSigma(Theta.tilde, tau.tilde, n)
    t2 <- Sys.time()
    print(difftime(t2,t1))
    
    S <- Sigma.list[[i]]
    Si <- S^(-1)
    
    Delta.list[[i]] <- new.Delta
    Delta <- new.Delta
  }
  
  new.Delta <- matrix(1, d, 1)
  Delta2 <- matrix(1,d,d)
  tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
  Theta.tilde <- computeTt(new.Delta, Theta.hat, m.lr.mat, lr.m.mat, m.shared.mat, l.ij.mat)
  Sigma.list[[d]] <- buildSigma(Theta.tilde, tau.tilde, n)
  
  Delta.list[[d]] <- new.Delta
  
  list(D = Delta.list,S = Sigma.list, Th = Tau.hat)
}