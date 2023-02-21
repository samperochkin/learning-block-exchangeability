sapply(list.files("~/Git/testing-tau/sim/sim-main/functionsLowBlock2/", full.names = T), source) 

image2 <- function(M, ...){
  image(t(M[ncol(M):1,]),...)
}

image3 <- function(M, ...){
  m <- round(c(M),5)
  vals <- sort(unique(m))
  
  for(k in seq_along(vals)){
    M[m == vals[k]] <- k
  }
  
  image(t(M[nrow(M):1,]), col=rainbow(length(vals)))
}

image4 <- function(M, d, only.k){
  
  if(only.k != F) M <- M[,1:only.k,drop=F]
  
  ij <- t(combn(d,2))
  par(mfrow=c(1,ncol(M)), mar = c(2,2,1,1))
  
  for(k in 1:ncol(M)){
    
    v <- M[,k]
    m <- round(c(v),5)
    vals <- sort(unique(m))
    
    V <- matrix(NA,d,d)
    V[ij] <- V[ij[,2:1]] <- m
    print(V)
    
    V[ij] <- V[ij[,2:1]] <- match(m,vals)
    
    image(t(V[d:1,]), col=rainbow(length(vals)))
  }
  
  par(mfrow=c(1,1))
}


# Setup -------------------------------------------------------------------

sig <- .7
n <- 200
clus <- c(1,1,1,1,2,2,2,2)
K <- length(unique(clus))

d <- length(clus)
p <- choose(d,2)

ij.mat <- t(combn(d,2))
kl.mat <- matrix(clus[ij.mat],p,2)
kl.mat2 <- unique(kl.mat)

# Compute Sigma
ts <- replicate(100, {
  X <- mvtnorm::rmvnorm(n,rep(0,d),(1-sig)+diag(d)*sig)
  pcaPP::cor.fk(X)[ij.mat]
})
S0 <- n*cov(t(ts))
S <- averageSigmaBlock(S0,ij.mat,clus)

# Eigen-decomposition
eigS <- eigen(S)
eigS$val
image(t(eigS$vec[p:1,]))

# B matrices
B0 <- matrix(0,d,K)
B0[cbind(1:d,clus)] <- 1
B <- apply(kl.mat2, 1, function(kl){
  tcrossprod(B0[,kl[1]],B0[,kl[2]])[ij.mat]
})
dim(B)
Bp <- MASS::ginv(B)

# eigen 1
R <- Bp %*% S %*% B
dim(R)
dim(S)

eigR <- eigen(R)
eigR$values
unique(round(eigS$val,5))

D <- diag(eigR$values)
V <- eigR$vec
R
V %*% D %*% solve(V)

W <- B %*% V
(S %*% W[,1])/W[,1]


# eigen 2
ee <- apply(ij.mat, 1, function(ij){
  clus[ij] <- 3
  K <- length(unique(clus))
  kl.mat <- matrix(clus[ij.mat],p,2)
  kl.mat2 <- unique(kl.mat)
  
  B0 <- matrix(0,d,K)
  B0[cbind(1:d,clus)] <- 1
  B <- apply(kl.mat2, 1, function(kl){
    tcrossprod(B0[,kl[1]],B0[,kl[2]])[ij.mat]
  })
  dim(B)
  Bp <- MASS::ginv(B)
  
  R <- Bp %*% S %*% B
  eigR <- eigen(R)
  round(eigR$val,5)
})

unique(round(eigS$values,5)) %in% unique(unlist(ee))

sort(unique(unlist(ee)), decreasing = T)
sort(unique(round(eigS$values,5)), decreasing = T)

table(unlist(ee))
table(round(eigS$values,5))





clus[c(1,8)] <- 3
K <- length(unique(clus))
kl.mat <- matrix(clus[ij.mat],p,2)
kl.mat2 <- unique(kl.mat)

B0 <- matrix(0,d,K)
B0[cbind(1:d,clus)] <- 1
B <- apply(kl.mat2, 1, function(kl){
  tcrossprod(B0[,kl[1]],B0[,kl[2]])[ij.mat]
})
dim(B)
Bp <- MASS::ginv(B)

R <- Bp %*% S %*% B
dim(R)
dim(S)

eigR <- eigen(R)
round(eigR$values,7)
unique(round(eigen(S)$val,7))

D <- diag(eigR$values)
V <- eigR$vec
R
V %*% D %*% solve(V)

W <- B %*% V
(S %*% W[,1])/W[,1]



table(round(eigen(S)$val,7))
