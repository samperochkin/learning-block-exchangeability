pathBuilder <- function(Tau.hat, Theta.hat, n, cutoff = NULL, w = .1){
  
  #### Package for parallel computations ####
  
  library(parallel)
  library(dendextend)
  
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
  tau.hat <- Tau.hat[l.ij.mat]
  
  #### Additional initial values ####
  
  Tau.list <- list()
  Sigma.list <- list()
  Delta.list <- list()
  
  Tau.list[[1]] <- Tau.hat
  tau.hat <- Tau.hat[l.ij.mat]
  
  Sigma.hat <- buildSigma(Theta.hat, tau.hat, n)
  Sigma.list[[1]] <- Sigma.hat
  
  S <- Sigma.list[[1]]
  S <- (1-w)*S + w*mean(S)
  Si <- S^(-1)
  
  Delta <- diag(d)
  Delta.list[[1]] <- Delta
  
  
  #### Functions for the iterations ####
  
  bigWinner <- function(Delta){
    K <- dim(Delta)[2]
    
    if(cutoff < d){
      Tau.tilde.K <- (t(Delta) %*% (Tau.hat - diag(d)) %*% Delta)/(t(Delta) %*% (1 - diag(d)) %*% Delta)
      hc <- hclust(as.dist(1 - abs(Tau.tilde.K)))
      
      anchors <- seq(1,K-1,ceiling((K-1)/20))
      cluster_mat <- sapply(anchors, function(k){
        cutree(hc, k)
      })
        
      total <- sapply(seq_along(anchors), function(k){
        sum(sapply(1:anchors[k], function(j){
          size <- sum(cluster_mat[,k] == j)
          size*(size-1)/2
        }))
      })

      cond <- total <= cutoff
      if(sum(cond) > 0){
        winner <- which.max(total*cond)
      }else{
        winner <- which.min(total)
      }
      
      cluster <- cluster_mat[,winner]

      ij.pairs <- do.call(rbind, 
                          lapply(1:anchors[winner], function(k){
                            members <- which(cluster == k)
                            if(length(members) == 1){
                              return(NULL)
                            }
                            t(combn(members,2))
                          }))
      L <- nrow(ij.pairs)

    }else{
      print("got there")
      L <- (K*(K-1)/2)
      ij.pairs <- t(sapply(1:L,function(l){
        i <- K-sum(l<=cumsum(K-(1:(K-1))))
        j <- l - ((i-1)*K - i*(i-1)/2 - i)
        c(i,j)
      }))
    }
  
    #cl <- makeCluster(detectCores()-1)
    #clusterExport(cl, c("Tau.hat", "tau.hat", "Si", "l.ij.mat", "d"), envir = environment())

    #the_gal <- which.min(parSapply(cl, 1:L, function(l){
    the_gal <- which.min(sapply(1:L, function(l){
      pair <- ij.pairs[l,]
      Delta[,pair[1]] <- rowSums(Delta[,pair])
      Delta <- Delta[,-pair[2]]
      Delta2 <- tcrossprod(Delta)
      tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
      vec <- tau.hat - tau.tilde
      crossprod(vec^2, Si)
    }))
    
    #stopCluster(cl)
    
    pair <- ij.pairs[the_gal,]
    Delta[,pair[1]] <- rowSums(Delta[,pair])
    Delta[,-pair[2]]
  }
  
  #### Iterations ####
  for(i in 2:(d-1)){
    t1.full <- Sys.time()

    #print(paste0("Hey! Now at iteration ", i))

    # We parallel the computations for the potential losses
    # cl <- makeCluster(detectCores()-1)
    # t1 <- Sys.time()
    # #print("export")
    # # Should definitively let the workers construct the candidates (they need) - to do later, or never.
    # # The bottleneck clearly this next line.
    # clusterExport(cl, c("Tau.hat", "tau.hat", "Si", "l.ij.mat", "d"), envir = environment())
    # t2 <- Sys.time()
    # #print(difftime(t2,t1))
    
    #print("comp")
    #t1 <- Sys.time()
    new.Delta <- bigWinner(Delta.list[[i-1]])
    #t2 <- Sys.time()
    #print(difftime(t2,t1))
    
    #stopCluster(cl)
    
    #print("housekeeping")
    #t1 <- Sys.time()
    Delta2 <- tcrossprod(new.Delta)
    tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
    Theta.tilde <- computeTt(new.Delta, Theta.hat, l.ij.mat)
    S <- buildSigma(Theta.tilde, tau.tilde, n)
    Sigma.list[[i]] <- S

    S <- (1-w)*S + w*mean(S)
    Si <- S^(-1)
    
    
    Delta.list[[i]] <- new.Delta
    Delta <- new.Delta
    
    #t2 <- Sys.time()
    #print(difftime(t2,t1))
    
    if(i %% 10 == 0){
      print(paste("full pass took",difftime(Sys.time(),t1.full), "--- Iteration", i))
    }
  }
  
  new.Delta <- matrix(1, d, 1)
  Delta2 <- matrix(1,d,d)
  tau.tilde <- ((Delta2 %*% (Tau.hat - diag(d)) %*% Delta2) / (Delta2 %*% (1 - diag(d)) %*% Delta2))[l.ij.mat]
  Theta.tilde <- computeTt(new.Delta, Theta.hat, l.ij.mat)
  Sigma.list[[d]] <- buildSigma(Theta.tilde, tau.tilde, n)
  
  Delta.list[[d]] <- new.Delta
  
  list(D = Delta.list,S = Sigma.list, Th = Tau.hat)
}
