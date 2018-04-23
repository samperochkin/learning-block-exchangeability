pathBuilder <- function(Tau.hat, Theta.hat, n, buffer = NULL, cutoff = 20, w = .1){
  
  #### Package for parallel computations ####
  #library(parallel)
  library(dendextend)
  
  #### fixed parameters ####
  d <- dim(Tau.hat)[1]
  p <- d * (d - 1) / 2
  
  #### vectorizing ####
  l.ij.mat <- t(sapply(1:p,function(l){
    i <- d-sum(l<=cumsum(d-(1:(d-1))))
    j <- l - ((i-1)*d - i*(i-1)/2 - i)
    c(i,j)
  }))
  
  ij.l.mat <- matrix(0,d,d)
  ij.l.mat[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(1:p,1:p)
  
  tau.hat <- Tau.hat[l.ij.mat]
  
  #### Additional initial values ####
  
  Tau.list <- list()
  Sigma.list <- list()
  clus.list <- list()
  
  Tau.list[[1]] <- Tau.hat
  Tt <- Tau.hat
  size.mat <- 1 - diag(d)
  tau.hat <- Tau.hat[l.ij.mat]
  tau.tilde <- tau.hat
  
  Sigma.hat <- buildSigma(Theta.hat, tau.hat, n)
  Sigma.list[[1]] <- Sigma.hat
  
  S <- Sigma.list[[1]]
  S <- (1-w)*S + w*mean(S)
  Si <- S^(-1)
  
  clus <- 1:d
  clus.list[[1]] <- clus
  
  #### Functions for the iterations ####
  
  reduceTt <- function(Tt,pair){
    if(nrow(Tt) <= 3){
      #print(nrow(Tt))
      Tt[-pair,pair[1]] <- sum(Tt[-pair,pair]*size.mat[-pair,pair])/sum(size.mat[-pair,pair])
    }else{
      Tt[-pair,pair[1]] <- rowSums(Tt[-pair,pair]*size.mat[-pair,pair])/rowSums(size.mat[-pair,pair])
    }
    Tt[pair,pair[1]] <- sum(Tt[pair,pair]*size.mat[pair,pair])/sum(size.mat[pair,pair])
    Tt[pair[1],] <- Tt[,pair[1]]
    as.matrix(Tt[-pair[2],-pair[2]])
  }
  reduceSM <- function(size.mat,pair){
    #print("in")
    if(nrow(size.mat) <= 3){
      size.mat[-pair,pair[1]] <- sum(size.mat[-pair,pair])
    }else{
      size.mat[-pair,pair[1]] <- rowSums(size.mat[-pair,pair])
    }
    size.mat[pair,pair[1]] <- sum(size.mat[pair,pair])
    size.mat[pair[1],] <- size.mat[,pair[1]]
    #print("out")
    size.mat[-pair[2],-pair[2]]
  }
  
  bigWinner <- function(clus){
    
    K <- length(unique(clus))
    
    if(!is.null(buffer)){
      hc <- hclust(as.dist(1 - abs(Tt)))
      
      #t1 <- Sys.time()
      anchors <- seq(max(1,K-buffer),K-1,ceiling(buffer/3))
      cluster_mat <- sapply(anchors, function(k){
        cutree(hc, k)
      })
      
      total <- sapply(seq_along(anchors), function(k){
        sum(sapply(1:anchors[k], function(j){
          size <- sum(cluster_mat[,k] == j)
          size*(size-1)/2
        }))
      })
      #t2 <- Sys.time()
      #print(paste("time1:",difftime(t2,t1)))
      
      
      cond <- total <= cutoff
      if(sum(cond) > 0){
        winner <- which.max(total*cond)
      }else{
        winner <- which.min(total)
      }
      
      #t1 <- Sys.time()
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
      #t2 <- Sys.time()
      #print(paste("time2:",difftime(t2,t1)))
      
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
    #t1 <- Sys.time()
    the_gal <- which.min(sapply(1:L, function(l){
      pair <- ij.pairs[l,]
      
      Tt <- reduceTt(Tt,pair)
      
      mem <- which(clus %in% pair)
      clus[mem] <- pair[1]
      clus <- frank(clus,ties.method="dense")
      #t1 <- Sys.time()
      ind <- do.call(what = rbind, lapply(seq_along(mem), function(i){
        cbind(rep(mem[i],d-i),(1:d)[-mem[1:i]])
      }))
      
      tau.tilde[ij.l.mat[ind]] <- Tt[cbind(clus[ind[,1]],clus[ind[,2]])]
      #t2 <- Sys.time()
      #print(paste("time3-inner:",difftime(t2,t1)))

      vec <- tau.hat - tau.tilde
      sum(vec^2*Si)
    }))
    
    #t2 <- Sys.time()
    #print(paste("time3:",difftime(t2,t1)))
    
    #stopCluster(cl)
    
    ij.pairs[the_gal,]
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
    winning.pair <- bigWinner(clus)
    mem <- which(clus %in% winning.pair)
    clus[mem] <- winning.pair[1]
    clus <- frank(clus,ties.method="dense")
    
    Tt <- reduceTt(Tt,winning.pair)
    size.mat <- reduceSM(size.mat,winning.pair)
    
    #t2 <- Sys.time()
    #print(difftime(t2,t1))
    
    #stopCluster(cl)
    
    #print("housekeeping")
    #t1 <- Sys.time()
    ind <- do.call(what = rbind, lapply(seq_along(mem), function(i){
      cbind(rep(mem[i],d-i),(1:d)[-mem[1:i]])
    }))
    tau.tilde[ij.l.mat[ind]] <- Tt[cbind(clus[ind[,1]],clus[ind[,2]])]
    #t2 <- Sys.time()
    #print(difftime(t2,t1))
    
    #Theta.tilde <- computeTt(Delta, Theta.hat, l.ij.mat)
    #S <- buildSigma(Theta.tilde, tau.tilde, n)
    Sigma.list[[i]] <- S
    
    S <- (1-w)*S + w*mean(S)
    Si <- S^(-1)
    
    clus.list[[i]] <- clus
    
    
    if(i %% 10 == 0){
      print(paste("full pass took",difftime(Sys.time(),t1.full), "--- Iteration", i, "of", d))
    }
  }
  
  #new.Delta <- matrix(1, d, 1)
  mem <- 1:d
  clus[mem] <- 1
  
  Tt <- reduceTt(Tt,c(1,2))
  size.mat <- reduceSM(size.mat,c(1,2))
  ind <- do.call(what = rbind, lapply(seq_along(mem), function(i){
    cbind(rep(mem[i],d-i),(1:d)[-mem[1:i]])
  }))
  tau.tilde[ij.l.mat[ind]] <- Tt[cbind(clus[ind[,1]],clus[ind[,2]])]
  
  #Theta.tilde <- computeTt(new.Delta, Theta.hat, l.ij.mat)
  #Sigma.list[[d]] <- buildSigma(Theta.tilde, tau.tilde, n)
  clus.list[[d]] <- clus
  
  list(clus = clus.list,S = Sigma.list, Th = Tau.hat)
}
