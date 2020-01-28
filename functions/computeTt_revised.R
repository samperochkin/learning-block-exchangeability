computeTt <- function(new.Delta, Theta.hat, m.lr.mat, lr.m.mat, m.shared.mat, l.ij.mat){
  
  d <- dim(new.Delta)[1]
  p <- d * (d-1) / 2
  K <- dim(new.Delta)[2]
  
  theta.hat <- Theta.hat[m.lr.mat]
  new.theta <- theta.hat
  
  
  clusters <- c(tcrossprod(1:K,new.Delta))
  clusters2 <- c(0,clusters)
  
  blocks2.mat <- matrix(clusters[c(l.ij.mat)], ncol = 2)
  blocks2.mat <- t(apply(FUN = sort, MARGIN = 1, X = blocks2.mat))
  b1 <- blocks2.mat[,1]
  blocks1.vec <- (b1 - 1) * (K - b1/2 + 1) + blocks2.mat[,2] - b1
  rm(b1)
  
  # create a list of blocks
  blocks.list <- list()
  set <- 1:p
  k <- 1
  while(length(set) > 0){
    ind <- which(blocks1.vec[set] == blocks1.vec[set[1]])
    blocks.list[[k]] <- set[ind]
    set <- set[-ind]
    k <- k + 1
  }
  
  bigger.than.1 <- which(sapply(blocks.list, length) > 1)
  
  K2 <- length(blocks.list)
  K3 <- length(bigger.than.1)
  
  # mapply on the list
  #(k1,k2) -> l1,l2,l3,l4  ,  r1,r2,r3,r4,r5,...
  # m.vec <- mapply(lr.m.mat[l,r])
  
  if(K == 1){
    block.pairs <- matrix(1,1,2)
  }else{
    
    within1 <- cbind(bigger.than.1, bigger.than.1)
    if(K3 > 1){
      between1 <- t(combn(bigger.than.1, 2))
    }else{
      between1 <- NULL
    }
    
    between2 <- cbind(
      rep(bigger.than.1,each = K2 - K3),
      rep((1:K2)[-bigger.than.1], times = K3)
    )
      
    block.pairs <- rbind( within1, between1, between2)
  }
  
  for(s in 1:dim(block.pairs)[1]){
    b <- block.pairs[s,]
    block1 <- blocks.list[[b[1]]]
    block2 <- blocks.list[[b[2]]]
    
    m.vec <- unique(c(lr.m.mat[block1, block2]))
    m.shared.groups <- matrix(clusters2[m.shared.mat[m.vec,] + 1], ncol = 2)
    m.shared.groups[m.shared.groups[,1] != 0] <- t(apply(FUN = sort, X = m.shared.groups[m.shared.groups[,1] != 0,], MARGIN = 1))
    
    set <- 1:length(m.vec)
    while(length(set) > 0){
      TF.vec <- sapply(set, function(m){
        identical(m.shared.groups[m,],m.shared.groups[set[1],])
      })
      new.theta[m.vec[set[TF.vec]]] <- mean(theta.hat[m.vec[set[TF.vec]]])
      set <- set[!TF.vec]
    }
  }
  
  new.Theta <- matrix(0, p, p)
  new.Theta[m.lr.mat] <- new.theta
  new.Theta <- new.Theta + t(new.Theta) - diag(diag(new.Theta)) 
  
  new.Theta
}

#image(t(new.Theta[p:1,]), zlim = c(min(new.Theta[p:1,]),max(new.Theta[p:1,])), axes=FALSE, col=colfunc(100))
#identical(new.Theta, Theta.list[[3]])
