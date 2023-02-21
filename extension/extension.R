library(data.table)
library(magrittr)

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
clus <- c(1,1,1,1,2,2,2,2,3,3,3,3)
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


table(round(eig$val,5))

#### inner block
d1 <- sum(clus==1)
B <- matrix(0,p,1)
B.star <- matrix(0,p,d)
id1 <- which(rowSums(kl.mat == 1) == 2)

B[id1] <- 1
B.star[cbind(id1,ij.mat[id1,1])] <- B.star[cbind(id1,ij.mat[id1,2])] <- 1
B.star <- B.star[,-which(colSums(B.star) == 0)]

G <- round(B %*% MASS::ginv(B),7)
G.star <- round(B.star %*% MASS::ginv(B.star),7)

I <- matrix(0,p,p)
I[cbind(id1,id1)] <- 1

IG.star <- I - G.star
GG.star <- G.star - G

V1 <- round(t(IG.star) %*% S,5)/IG.star
V2 <- round(t(GG.star) %*% S,5)/GG.star

unique(round(V1[!is.nan(V1)],5))
unique(round(V2[!is.nan(V2) & !is.infinite(V2)],5))
table(round(eigen(S)$val,6))

v <- IG.star[,1]
unique(round(S %*% v,5)/v)

v <- GG.star[,1]
unique(round(S %*% v,5)/v)



V1 <- matrix(0,d,d)
V2 <- matrix(0,d,d)
V1[ij.mat] <- V1[ij.mat[,2:1]] <- IG.star[,2]
V2[ij.mat] <- V2[ij.mat[,2:1]] <- GG.star[,2]
image(V1)
image(V2)

ks <- which(colSums(abs(GG.star)) != 0)
v <- rowSums(GG.star)

fun <- function(as){
  V <- matrix(0,d,d)
  V[1,i1] <- as[1]
  V[2,i1] <- as[1]
  V[1,i2] <- as[2]
  V[2,i2] <- as[2]
  V[1,i3] <- as[3]
  V[2,i3] <- as[3]
  sd((S[,1:(2*d-3)] %*% V[ij.mat[1:(2*d-3),]]) - as[5]*V[ij.mat[1:(d-1),]]) + 1000*abs(1-crossprod(V[ij.mat[1:(d-1),]])) + 1000*abs(sum(V))
}
fun <- function(as){
  V <- matrix(0,d,d)
  V[1,i1] <- as[1]
  V[1,i2] <- as[2]
  V[1,i3] <- as[3]
  sd((S[,1:(d-1)] %*% V[ij.mat[1:(d-1),]]) - as[4]*V[ij.mat[1:(d-1),]]) + 1000*abs(1-crossprod(V[ij.mat[1:(d-1),]])) + 1000*abs(sum(V))
}

opt <- optim(par = rep(1,5), fun)
as <- opt$par
V <- matrix(0,d,d)
V[1,i1] <- as[1]
V[1,i2] <- as[2]
V[1,i3] <- as[3]
V[2,i1] <- as[4]
V[2,i2] <- as[5]
V[2,i3] <- as[6]
crossprod(V)
(S %*% V[ij.mat]) / V[ij.mat]

##### outer block
d1 <- sum(clus==1)
d2 <- sum(clus==2)
B <- matrix(0,p,1)
B.star <- matrix(0,p,d)
id12 <- which(kl.mat[,1] == 1 & kl.mat[,2] == 2)

B[id12] <- 1
B.star[cbind(id12,ij.mat[id12,1])] <- B.star[cbind(id12,ij.mat[id12,2])] <- 1
B.star <- B.star[,-which(colSums(B.star) == 0)]

G <- round(B %*% MASS::ginv(B),8)
G.star <- round(B.star %*% MASS::ginv(B.star),8)

I <- matrix(0,p,p)
I[cbind(id12,id12)] <- 1

IG.star <- I - G.star
GG.star <- G.star - G

(round(t(IG.star) %*% S,8)/IG.star)[1:10,1:10]
(round(t(GG.star) %*% S,8)/GG.star)[1:10,1:10]
table(round(eigen(S)$val,5))

v <- IG.star[,6]
unique(round(S %*% v,8)/v)

v <- GG.star[,4]
unique(round(S %*% v,8)/v)

V1 <- matrix(0,d,d)
V2 <- matrix(0,d,d)
V1[ij.mat] <- V1[ij.mat[,2:1]] <- IG.star[,4]
V2[ij.mat] <- V2[ij.mat[,2:1]] <- GG.star[,4]
image(V1)
image(V2)






#### inner block
i1 <- which(clus == 1)
d1 <- length(i1)
i2 <- which(clus == 2)
i3 <- which(clus == 3)
B <- matrix(0,p,3)
B.star <- matrix(0,p,3)

id1 <- which(rowSums(kl.mat == 1) == 2)
id2 <- which(kl.mat[,1] == 1 & kl.mat[,2] == 2)
id3 <- which(kl.mat[,1] == 1 & kl.mat[,2] == 3)

sapply(id1, function(r){
  ij <- ij.mat[r,]
  
})


B[id1,1] <- 1
B[id2,2] <- 1
B[id3,3] <- 1

B.star1[cbind(id1,ij.mat[id1,1])] <- B.star1[cbind(id1,ij.mat[id1,2])] <- 1
B.star2[cbind(id2,ij.mat[id2,1])] <- B.star2[cbind(id2,ij.mat[id2,2])] <- 1
B.star3[cbind(id3,ij.mat[id3,1])] <- B.star3[cbind(id3,ij.mat[id3,2])] <- 1


G <- round(B %*% MASS::ginv(B),7)
G.star1 <- round(B.star1 %*% MASS::ginv(B.star1),7)
G.star2 <- round(B.star2 %*% MASS::ginv(B.star2),7)
G.star3 <- round(B.star3 %*% MASS::ginv(B.star3),7)
G.star <- G.star1 + G.star2 + G.star3

GG.star <- G.star - G

V2 <- round(t(GG.star) %*% S,5)/GG.star

unique(round(V1[!is.nan(V1)],5))
unique(round(V2[!is.nan(V2) & !is.infinite(V2)],5))
table(round(eigen(S)$val,6))

v <- IG.star[,1]
unique(round(S %*% v,5)/v)

v <- GG.star[,1]
unique(round(S %*% v,5)/v)
