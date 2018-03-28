
delta <- (Delta %*% (1 - diag(6)) %*% ginv(Delta))[l.ij.mat]
cbind((ddd %*% t(ddd) - diag(p)) %*% tau.hat, tau.hat)

DD <- tcrossprod(Delta)

delta <- (DD/(DD %*% (1 - diag(d)) %*% DD))[l.ij.mat]


tt <- (DD %*% (Tau.hat - diag(d)) %*% DD / (DD %*% (1 - diag(d)) %*% DD))[l.ij.mat]
tau <- Tau[l.ij.mat]

cbind(tau,tau.hat,tt)

cbind(tau, tt, delta %*% t(delta) %*% tau.hat)


delta1 <- tcrossprod(Delta)[l.ij.mat]
delta1 <- ginv(tcrossprod(Delta))[l.ij.mat]
delta2 <- (2/(DD %*% (1 - diag(d)) %*% DD))[l.ij.mat]

delta1 <- (1/((1 - diag(d)) %*% DD))[l.ij.mat]
delta2 <- ginv(DD)[l.ij.mat]

cbind((delta1 %*% t(delta2)) %*% tau.hat, tt)

mat <- matrix(colSums(Delta), nrow = 6, ncol = 6)
mat <- (mat * t(mat) - diag(colSums(Delta)))
diag(mat) <- diag(mat)/2

delta <- (Delta %*% mat %*% t(Delta))[l.ij.mat]

delta1 <- tcrossprod(Delta)[l.ij.mat]
delta2 <- tcrossprod(Delta)[l.ij.mat[,2:1]]

delta <- (tcrossprod(Delta) %*% ginv(tcrossprod(Delta)))[l.ij.mat]

matrix(delta, p, p)

(delta %*% t(delta))[1:10,1:10]
