computeThPara <- function(X){
  
  d <- ncol(X)
  p <- d * (d-1) / 2
  n <- nrow(X)
  
  l.ij.mat <- t(combn(1:d,2))

  indicatorMatrixFun <- function(Y){
    M1 <- matrix(Y, n, n, byrow = TRUE)
    M2 <- matrix(Y, n, n)
    1 * (M1 < M2)
  }
  
  l.unique <- 1:p
  i.needed <- 1:d

  print("Here?")
  
  I.list <- sapply(1:d,function(i){indicatorMatrixFun(X[,i])},simplify=FALSE)
  #I2.list <- sapply(1:p,function(l){i <- l.ij.mat[l,1]; j <- l.ij.mat[l,2]; I.list[[i]]*I.list[[j]]}, simplify = FALSE)

  #IJ.list <- sapply(1:p, function(ind){I <- I2.list[[ind]]; I + t(I)}, simplify = FALSE)
  one.vec <- rep(1,n)
  
  print("or there?")

  #IJ1.vec <- sapply(1:p, function(l){t(one.vec) %*% IJ.list[[l]]})
  #IJ2.vec <- sapply(1:p, function(l){IJ.list[[l]] %*% one.vec})

  sapply(1:p, function(r){
    
    i <- l.ij.mat[r,1]
    j <- l.ij.mat[r,2]
    I2 <- I.list[[i]]*I.list[[j]]
    IJ <- I2 + t(I2)
    
    if(r %% 5000 == 0){
      print(c("progress:",r/p))
    }
    #(4 / (n * (n - 1)))^2 * (tcrossprod((t(one.vec) %*% IJ.list[[r]])) - crossprod(c(I2.list[[r]]),c(IJ.list[[r]])))
    (4 / (n * (n - 1)))^2 * (tcrossprod(t(one.vec) %*% IJ) - crossprod(c(I2),c(IJ)))
  })
}
