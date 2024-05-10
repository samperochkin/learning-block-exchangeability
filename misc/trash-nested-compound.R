


# packages and functions --------------------------------------------------

# install_github("samperochkin/tautest")
library(tautests)


constrainSigma <- function(Sh, clus, ij = NULL, debug = F){
  
  #### fixed parameters ####
  d <- length(clus)
  p <- d * (d - 1) / 2
  K <- length(unique(clus))
  
  #### vectorizing ####
  
  if(is.null(ij)) ij <- tautests::R2IJ(1:p)
  kl <- matrix(clus[c(ij)],p,2) |> apply(1, sort) |> t()
  kl0 <- unique(kl)
  L <- nrow(kl0)
  
  s_kl <- kl[,1] + (K+1)*kl[,2]
  s_kl0 <- unique(s_kl) # should be in same order as u_kl
  s_kl <- match(s_kl, s_kl0)
  if(length(s_kl0) != nrow(kl0)) cat("problem in generating unique block ids\n")
  
  kl0_type <- apply(kl0, 1, function(x) length(unique(x)))

  rr <- expand.grid(1:L, 1:L)
  klkl0_type <- cbind(kl0_type[rr[,1]],kl0_type[rr[,2]],
                      apply(rr, 1, function(x) length(intersect(kl0[x[1],],kl0[x[2],]))))
  klkl0_type <- t(apply(klkl0_type, 1, sort))
  klkl0_type_s <- 2^klkl0_type[,1] * 3^klkl0_type[,2] * 5^klkl0_type[,3]
  klkl0_type_s <- match(klkl0_type_s, unique(klkl0_type_s))
  
  # loop to make sure we don't run into memory problems
  for(m in unique(klkl0_type_s)){
    
    # find all pairs of block of this type
    rr0 <- rr[klkl0_type_s == m,,drop=F]
    
    rs <- do.call("rbind", lapply(1:nrow(rr0), function(q){
      r_u <- which(s_kl == rr0[q,1])
      s_v <- which(s_kl == rr0[q,2])
      
      ij_u <- ij[r_u,,drop=F]
      ij_v <- ij[s_v,,drop=F]
      p_u <- nrow(ij_u)
      p_v <- nrow(ij_v)
      
      psi <- mapply(FUN = function(x,y){
        psi <- length(intersect(ij_u[x,],ij_v[y,]))
        c(rep(0,2-length(psi)), sort(psi))
      }, x=rep(1:p_u, times = p_v), y=rep(1:p_v, each = p_u)) |>
        t() |> c() |> array(dim = c(p_u, p_v, 2))
      
      s_psi <- psi[,,1] + (L+1)*psi[,,2]
      groups <- split(1:length(s_psi), s_psi)
      
      lapply(groups, function(g){
        s <- ceiling(g/p_u)
        r <- g-(s-1)*p_u
        cbind(r_u[r],s_v[s])
      })
    }))
    
    for(w in 1:ncol(rs)){
      rs0 <- do.call("rbind", rs[,w])
      print(cbind(ij[rs0[,1],],ij[rs0[,2],]))
      Sh[rs0] <- Sh[rs0[,2:1]] <- mean(Sh[rs0])
    }
  }

  return(Sh)
}


# data and basic quantities
n <- 100
d <- 12
p <- choose(d,2)
clus <- rep(1:2, each = d/2)
# clus <- c(1,1,1,1,2,2,2,2,3,3)
K <- length(unique(clus))

# ij <- tautests::R2IJ(1:p)
ij <- do.call("rbind", c(
  lapply(1:K, function(k) t(combn((k-1)*d/K + 1:(d/K),2))),
  lapply(1:(K-1), function(k){
    a1 <- (k-1)*d/K; a2 <- 1:(d/K)
    as.matrix(expand.grid(a1+a2, (a1+d/K+1):d))
  })
))
  
  

sigs <- runif(K,.1,.25)
sigs <- rep(.1,K)
Z <- sapply(sigs, function(s) rnorm(n,0,s))
z <- rnorm(n,0,.1)
X <- matrix(rnorm(n*d,0,.1),n,d) + Z[,clus] + z
TSh <- tautests::tau_and_jack(X, ijs = ij)
Sh <- TSh$jack_var
Th <- TSh$tau
th <- Th[ij]

St <- constrainSigma(Sh, clus, ij)
image2(St)
eig <- eigen(St)
image2(eig$vectors)
# t(eig$vectors) %*% St %*% eig$vectors
plot(eig$values); min(eig$values); table(round(eig$values, 5))
Sti <- eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
Sti2 <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)


A <- matrix(0,d,K)
A[cbind(1:d,clus)] <- 1  
A <- t(apply(ij, 1, function(x) colSums(A[x,])))
B <- matrix(0, nrow(A), 2)
B[cbind(1:p,apply(A,1,max))] <- 1
L <- ncol(B)


G <- B %*% MASS::ginv(B)
G <- (G + t(G))/2
image2(G)
# IG <- diag(p) - G
tt <- G %*% th
Tt <- diag(d); Tt[ij] <- Tt[ij[,2:1]] <- tt
par(mfrow=c(1,2), mar=c(2,2,1,1))
image(t(Th)[,d:1], zlim = range(c(Th,Tt)))
image(t(Tt)[,d:1], zlim = range(c(Th,Tt)))

b2 <- matrix(0, p, d); b2[cbind(1:p,ij[,1])] <- b2[cbind(1:p,ij[,2])] <- 1
A <- matrix(0,d,K)
A[cbind(1:d,clus)] <- 1  
A <- t(apply(ij, 1, function(x) colSums(A[x,])))
A_u <- unique(A)
C <- matrix(0, nrow(A), nrow(A_u))
C[cbind(1:p, apply(A, 1, function(a) which.max(colSums(a == t(A_u)))))] <- 1
C <- C[, colSums(C) != 0]
B2 <- cbind(B, C, b2)
G2 <- B2 %*% MASS::ginv(B2)
G2 <- (G2 + t(G2))/2
image2(G)
image2(G2)
matrixcalc::is.symmetric.matrix(G)
matrixcalc::is.symmetric.matrix(G2)

range((G2 %*% G) - G)
th[1:5]
(G %*% th)[1:5]
(G2 %*% th)[1:5]

ev <- eigen(St)$val
table(round(ev, 5))

G1 <- G[,1]
SG <- (St %*% G1)
xx <- SG[abs(G1) > 1e-10]/G1[abs(G1) > 1e-10]
xx
yy <- sapply(xx, function(x) ev[which.min(abs(ev - x))])
yy
cbind(xx, yy, xx - yy)


GG1 <- (G2 - G)[,1]
SG <- (St %*% GG1)
xx <- SG[abs(GG1) > 1e-10]/GG1[abs(GG1) > 1e-10]
xx
yy <- sapply(xx, function(x) ev[which.min(abs(ev - x))])
yy
plot(xx, yy)
abline(a = 0, b=1)

GG1[1:p <= 31] <- 0


IG1 <- (diag(p) - G2)[,1]
SG <- (St %*% IG1)
xx <- SG[abs(IG1) > 1e-10]/IG1[abs(IG1) > 1e-10]
sapply(xx, function(x) ev[which.min(abs(ev - x))])
plot(xx, sapply(xx, function(x) ev[which.min(abs(ev - x))]))
abline(a = 0, b=1)
image2(G2)
image2(St)


H <- B %*% solve(t(B) %*% Sti %*% B) %*% t(B) %*% Sti
image2(H)
range(H - G)

H2 <- B2 %*% MASS::ginv(t(B2) %*% Sti %*% B2) %*% t(B2) %*% Sti
image2(G2)
image2(H2)
range(H2 - G2)
range(H2 - t(H2))
image2(H2 - t(H2))
table(round(H2 - t(H2), 5))
image2(H2 - G2)



# test with eigen stuff
image2(St)
eig <- eigen(St)
ev <- round(eig$values, 5)
plot(ev); min(ev); table(ev)
ev_u <- unique(ev)
Sti <- eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
Sti2 <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)


G3 <- eig$vectors*matrix(ev == ev_u[9],p,p,byrow = T)
image2(G3)

image2(G)
image2(G2)
image2(G3)
matrixcalc::is.symmetric.matrix(G)
matrixcalc::is.symmetric.matrix(G2)

range((G3 %*% G) - G)
th[1]
(G %*% th)[1]
(G2 %*% th)[1]
(G3 %*% th)[1]

k <- which(colSums(G3) != 0)[1]
GG1 <- (G - G3)[,k]
SG <- (St %*% GG1)
xx <- SG[abs(GG1) > 1e-10]/GG1[abs(GG1) > 1e-10]
xx
yy <- sapply(xx, function(x) ev[which.min(abs(ev - x))])
yy





par(mfrow = c(2,5))
for(k in 1:length(ev_u)){
  G3 <- eig$vectors*matrix(ev == ev_u[k],p,p,byrow = T)
  ttt <- unique(round(Sti2 %*% (sqrt(n)*G3 %*% th), 7))
  hist(ttt[abs(ttt) > 1e-7], breaks = 20)
}


ev <- round(eig$values, 5)
ev_u <- unique(ev)
cols <- colorRampPalette(c("blue", "white", "red"))(200)
G3 <- eig$vectors[,ev == ev_u[3], drop=F]
nc <- ncol(G3)
zz <- max(abs(G3))*c(-1,1)
par(mfrow=c(1,nc+1))
image(t(G3[p:1,]), zlim = zz, col=cols)
for(k in 1:nc) image(t(G3[p:1,k,drop=F]), zlim = zz, col=cols)


zz <- max(abs(eig$vectors))*c(-1,1)
G3 <- eig$vectors[,sapply(ev_u, function(ee) which(ev == ee)[1])]
nc <- ncol(G3)
par(mfrow=c(1,nc))
for(k in 1:nc) image(t(G3[p:1,k,drop=F]), zlim = zz, col=cols)
