pathBuilder0Para <- function(Tau.hat, Theta.hat, n, cutoff = NULL){
  
  #### Package for parallel computations ####
  
  library(parallel)
  
  #### fixed parameters ####
  
  d <- dim(Tau.hat)[1]
  p <- d * (d - 1) / 2
  
  if(is.null(cutoff)){
    cutoff <- d
  }
  
  #### vectorizing ####
  
  l.ij.mat <- t(sapply(1:p,function(l){
    i <- d-sum(l<=cumsum(d-(1:(d-1))))
    j <- l - ((i-1)*d - i*(i-1)/2 - i)
    c(i,j)
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
    
    if(cutoff < d){
      Tau.tilde.K <- (t(Delta) %*% (Tau.hat - diag(d)) %*% Delta)/(t(Delta) %*% (1 - diag(d)) %*% Delta)
      hc <- hclust(as.dist(1 - abs(Tau.tilde.K)))
      
      cluster_mat <- sapply(1:K, function(j){
        dendextend::cutree(hc, j)
      })
      
      clusters <- cluster_mat[,which(sapply(1:K, function(k){
        sum(sapply(1:k, function(j){
          size <- sum(cluster_mat[,k] == j)
          size*(size-1)/2
        }))
      }) <= cutoff)[1]]
      
      ij.pairs <- do.call(rbind, 
                          sapply(1:length(unique(clusters)), function(k){
                            if(length(which(clusters == k)) == 1){
                              return(NULL)
                            }
                            t(combn(which(clusters == k),2))
                          }, simplify = FALSE))
      L <- nrow(ij.pairs)
      
    }else{
      L <- (K*(K-1)/2)
      ij.pairs <- t(sapply(1:L,function(l){
        i <- K-sum(l<=cumsum(K-(1:(K-1))))
        j <- l - ((i-1)*K - i*(i-1)/2 - i)
        c(i,j)
      }))
    }
    
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
    vec <- Tau.hat[l.ij.mat] - tau.tilde
    crossprod(vec^2, Si)
  }
  
  #### Iterations ####
  # i <- 3
  
  for(i in 2:(d-1)){
    print(paste0("Hey! Now at iteration ", i))
    Delta.candidates <- DeltaCandidates(Delta.list[[i-1]])
    
    # We parallel the computations for the potential losses
    cl <- makeCluster(detectCores()-1)
    t1 <- Sys.time()
    print("export")
    # Should definitively let the workers construct the candidates (they need) - to do later, or never.
    # The bottleneck clearly this next line.
    clusterExport(cl, c("Delta.candidates", "Tau.hat", "Si", "l.ij.mat", "d", "loss", "tau.hat"), envir = environment())
    t2 <- Sys.time()
    print(difftime(t2,t1))
    
    new.Delta <- Delta.candidates[[which.min(parSapply(cl, Delta.candidates, function(Delta){loss(Delta, Tau.hat, Si, l.ij.mat, d)}))]]
    
    stopCluster(cl)
    
    print("rest")
    Delta2 <- tcrossprod(new.Delta)
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    Theta.tilde <- computeTt(new.Delta, Theta.hat, l.ij.mat)
    Sigma.list[[i]] <- buildSigma(Theta.tilde, tau.tilde, n)
    
    S <- Sigma.list[[i]]
    Si <- S^(-1)
    
    Delta.list[[i]] <- new.Delta
    Delta <- new.Delta
  }
  
  new.Delta <- matrix(1, d, 1)
  Delta2 <- matrix(1,d,d)
  tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
  Theta.tilde <- computeTt(new.Delta, Theta.hat, l.ij.mat)
  Sigma.list[[d]] <- buildSigma(Theta.tilde, tau.tilde, n)
  
  Delta.list[[d]] <- new.Delta
  
  list(D = Delta.list,S = Sigma.list, Th = Tau.hat)
}
