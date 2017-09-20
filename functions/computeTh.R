computeTh <- function(X){
  d <- dim(X)[2]
  p <- d * (d - 1) / 2
  n <- dim(X)[1]
  
  l.ij.mat <- t(sapply(1:p,function(l){
    i <- d-sum(l<=cumsum(d-(1:(d-1))))
    j <- l - ((i-1)*d - i*(i-1)/2 - i)
    c(i,j)
  }))
  
  indicatorMatrixFun <- function(Y){
    M1 <- matrix(Y, n, n, byrow = TRUE)
    M2 <- matrix(Y, n, n)
    1 * (M1 < M2)
  }
  
  I.list <- sapply(1:d,function(p){indicatorMatrixFun(X[,p])},simplify=FALSE)
  I2.list <- sapply(1:p,function(l){i <- l.ij.mat[l,1]; j <- l.ij.mat[l,2]; I.list[[i]]*I.list[[j]]}, simplify = FALSE)
  
  IJ.list <- sapply(1:p, function(l){I <- I2.list[[l]]; I + t(I)}, simplify = FALSE)
  one.vec <- rep(1,n)
  IJ1.vec <- sapply(1:p, function(l){t(one.vec) %*% IJ.list[[l]]})
  IJ2.vec <- sapply(1:p, function(l){IJ.list[[l]] %*% one.vec})
  
  I2.mat <- matrix(unlist(I2.list), n^2, p)
  IJ.mat <- matrix(unlist(IJ.list), n^2, p)
  
  (4 / (n * (n - 1)))^2 * (crossprod(IJ1.vec,IJ2.vec) - crossprod(I2.mat,IJ.mat))
}