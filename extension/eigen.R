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
clus <- c(1,1,1,1,1,2,2,2,2,3,3)
d <- length(clus)
p <- choose(d,2)
ij.mat <- t(combn(d,2))
kl.mat <- matrix(clus[ij.mat],p,2)

# Compute Sigma
ts <- replicate(100, {
  X <- mvtnorm::rmvnorm(n,rep(0,d),(1-sig)+diag(d)*sig)
  pcaPP::cor.fk(X)[ij.mat]
})
S0 <- n*cov(t(ts))
S <- averageSigmaBlock(S0,ij.mat,clus)

# Eigen-decomposition
eig <- eigen(S)
eig$val
image(t(eigen(S)$vec[p:1,]))

cbind(round(eig$vectors,5),ij.mat,kl.mat)
cbind(round(eig$vectors[,1],4),ij.mat,kl.mat)

image(abs(round(t(eigen(S)$vec[p:1,]),6)) < .000001)


# Analysis of non-zero entries --------------------------------------------

# encode groups (k,l) from which non-zero entries (i,j) come from
pat <- lapply(1:p, function(r){
  unique(kl.mat[abs(eig$vec[,r]) > .00001,])
})

invisible(
sapply(unique(pat), function(pa){
  
  print(pa)
  print(sum(sapply(pat, function(pa2) identical(pa,pa2))))
  print(table(round(eig$values[sapply(pat, function(pa2) identical(pa,pa2))],5)))
  print("--------------------")
}))

clus <- match(pat, unique(pat))

# K(K+1)/2 eigenvectors with no non-zero entries

sapply(1:max(clus), function(x) image4(eig$vectors[,clus == x], d, only.k = 1))
#
image4(S[,5,drop=F], d, only.k = 1)
#
image4(eig$vectors[,clus == 1], d, only.k = F)
image4(eig$vectors[,clus == 3], d, only.k = 3)
#






# v1
v1 <- function(a){
  a1 <- a[1]
  a2 <- a[2]
  a3 <- a[3]
  sapply(1:p, function(r) a1*identical(kl.mat[r,],c(1,1)) +
           a2*identical(kl.mat[r,],c(2,2)) +
           a3*identical(kl.mat[r,],c(1,2)))
}
fun <- function(a){
  sd((S %*% v1(a))/v1(a))
}

a <- optim(par = c(1,2,3), fn = fun)$par
v1(a)
lam <- (S %*% v1(a))/v1(a)
lam


# v2
s <- 1
v2 <- function(a){
  
  I <- function(r){ identical(kl.mat[r,],c(1,1)) * c(all(ij.mat[r,] == c(1,2)),
                                                     all(ij.mat[r,] == c(1,3)),
                                                     all(ij.mat[r,] == c(1,4)),
                                                     all(ij.mat[r,] == c(2,3)),
                                                     all(ij.mat[r,] == c(2,4)),
                                                     all(ij.mat[r,] == c(3,4))) }

    # identical(kl.mat[r,],c(2,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),ij.mat[r,]),
    # identical(kl.mat[r,],c(2,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),ij.mat[r,1]),
    # identical(kl.mat[r,],c(2,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),ij.mat[r,2]),
    # identical(kl.mat[r,],c(2,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),integer(0)),
    # identical(kl.mat[r,],c(1,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),ij.mat[r,]),
    # identical(kl.mat[r,],c(1,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),ij.mat[r,1]),
    # identical(kl.mat[r,],c(1,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),ij.mat[r,2]),
    # identical(kl.mat[r,],c(1,2))*identical(intersect(ij.mat[r,],ij.mat[s,]),integer(0)))}
  sapply(1:p, function(r) crossprod(a, I(r)))
}
fun <- function(a){
  v <- v2(a)
  res <- (S %*% v)/v
  sd(res[v != 0]) + abs(c(crossprod(v2(a))) - 1)
}

a <- optim(par = rep(1,6), fn = fun)$par
v2(a)
lam <- (S %*% v2(a))/v2(a)
lam




# v2
v2 <- function(a){
  a1 <- a[1]
  a2 <- a[2]
  a3 <- a[3]
  a4 <- a[4]
  sapply(1:p, function(r){
    a1*identical(kl.mat[r,],c(1,1)) +
      a2*identical(kl.mat[r,],c(2,2)) +
      a3*identical(kl.mat[r,],c(1,2))*(1 %in% ij.mat[r,]) +
      a4*(identical(kl.mat[r,],c(1,2)) & !(1 %in% ij.mat[r,]))})
}
fun <- function(a){
  sd((S %*% v2(a))/v2(a)) #+ abs(c(crossprod(v2(a)))-1)
}

a <- optim(par = c(1,2,3,4), fn = fun)$par
v2(a)
lam <- (S %*% v2(a))/v2(a)
lam




block <- match(kl.mat[,1] + kl.mat[,2]*100, unique(kl.mat[,1] + kl.mat[,2]*100))
B <- matrix(0,p,3)
B[cbind(1:p,block)] <- 1
BB <- B %*% MASS::ginv(B)

cbind(block,BB[,1])






S1 <- averageSigmaBlock(S0, ij.mat, clus = rep(1,d))

B1 <- rep(1,p)
# B1.star <- sapply(1:d, function(i) as.integer(rowSums(ij.mat == i) > 0))
B1.star <- sapply(1:d, function(i) as.integer(rowSums(ij.mat == i) > 0))

V1 <- G1 <- B1 %*% MASS::ginv(B1)
G1.star <- B1.star %*% MASS::ginv(B1.star)
V2 <- G1.star - G1

(V1 %*% S1)/V1
(V2 %*% S1)/V2



S2 <- averageSigmaBlock(S0, ij.mat, clus = clus)

B2 <- matrix(0,p,3)
B2[cbind(1:p,block)] <- 1
# B2.star <- sapply(1:d, function(i) as.integer(rowSums(ij.mat == i) > 0))
# B2.star <- cbind(sapply(1:d, function(i) as.integer(rowSums(ij.mat == i) > 0)),sapply(1:K, function(k) as.integer(rowSums(kl.mat == k) > 0)))
B2.star <- cbind(sapply(1:d, function(i) as.integer(rowSums(ij.mat == i) > 0)),sapply(1:K, function(k) as.integer(rowSums(kl.mat == k) > 0)))

V1 <- G2 <- B2 %*% MASS::ginv(B2)
G2.star <- (B2.star %*% MASS::ginv(B2.star))
V2 <- G2.star - G2

(V1 %*% S2)/V1
(V2 %*% S2)/V2



# S2 <- averageSigmaBlock(S0, ij.mat, clus = clus)
# 
# B2 <- sapply(1:K, function(k) as.integer(rowSums(kl.mat == k) > 0))
# V1 <- G2 <- B2 %*% MASS::ginv(B2)
# V1 <- apply(V1,1,function(v) v/c(crossprod(v)))
# (V1 %*% S2)/V1


