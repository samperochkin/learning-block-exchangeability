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
  
  I.list <- sapply(1:d,function(i){indicatorMatrixFun(X[,i])},simplify=FALSE)
  I2.list <- sapply(1:p,function(l){i <- l.ij.mat[l,1]; j <- l.ij.mat[l,2]; I.list[[i]]*I.list[[j]]}, simplify = FALSE)

  IJ.list <- sapply(1:p, function(ind){I <- I2.list[[ind]]; I + t(I)}, simplify = FALSE)
  one.vec <- rep(1,n)

  IJ1.vec <- sapply(1:p, function(l){t(one.vec) %*% IJ.list[[l]]})
  IJ2.vec <- sapply(1:p, function(l){IJ.list[[l]] %*% one.vec})

  I2.mat <- matrix(unlist(I2.list), n^2, p)
  IJ.mat <- matrix(unlist(IJ.list), n^2, p)

  sapply(1:p, function(r){
    (4 / (n * (n - 1)))^2 * (crossprod(IJ1.vec[,r],IJ2.vec[,r]) - crossprod(I2.mat[,r],IJ.mat[,r]))
  })
}