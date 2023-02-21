


# packages and functions --------------------------------------------------

# install_github("samperochkin/tautest")
library(tautests)

constrainSigma <- function(Sh, clus, debug = F){
  
  #### fixed parameters ####
  d <- length(clus)
  p <- d * (d - 1) / 2
  K <- length(unique(clus))
  
  #### vectorizing ####
  
  ij <- tautests::R2IJ(1:p)
  kl <- matrix(clus[c(ij)],p,2) |> apply(1, sort) |> t()
  kl0 <- unique(kl)
  L <- nrow(kl0)
  
  s_kl <- kl[,1] + (K+1)*kl[,2]
  s_kl0 <- unique(s_kl) # should be in same order as u_kl
  s_kl <- match(s_kl, s_kl0)
  if(length(s_kl0) != nrow(kl0)) cat("problem in generating unique block ids\n")
  
  # loop to make sure we don't run into memory problems
  for(u in 1:L){
    for(v in u:L){
      
      kl_u <- kl0[u,]
      kl_v <- kl0[v,]
      
      if(debug) cat("----------\n case", kl_u, "and", kl_v, "\n")
      
      r_u <- which(s_kl == u)
      s_v <- which(s_kl == v)
      
      if(length(intersect(kl_u,kl_v)) == 0){
        # no overlap possible
        Sh[r_u,s_v] <- Sh[s_v,r_u] <- mean(Sh[r_u,s_v])
        next
      }
      
      # collect entries from both blocks
      ij_u <- ij[r_u,]
      ij_v <- ij[s_v,]
      p_u <- nrow(ij_u)
      p_v <- nrow(ij_v)
      
      psi <- mapply(FUN = function(x,y){
        psi <- clus[intersect(ij_u[x,],ij_v[y,])]
        c(rep(0,2-length(psi)), sort(psi))
      }, x=rep(1:p_u, times = p_v), y=rep(1:p_v, each = p_u)) |>
        t() |> c() |> array(dim = c(p_u, p_v, 2))
      
      s_psi <- psi[,,1] + (L+1)*psi[,,2]
      groups <- split(1:length(s_psi), s_psi)
      
      for(g in groups){
        s <- ceiling(g/p_u)
        r <- g-(s-1)*p_u
        Sh[cbind(r_u[r],s_v[s])] <- Sh[cbind(s_v[s],r_u[r])] <- mean(Sh[cbind(r_u[r],s_v[s])])
        # image(t(Sh)[rev(r_u),s_v])
      }
      
      if(debug){
        l1 <- min(length(unique(kl_u)),length(unique(kl_v)))
        l2 <- max(length(unique(kl_u)),length(unique(kl_v)))
        l12 <- length(intersect(kl_u,kl_v))
        
        category <- which(c(l1 == 1 & l2 == 1 & l12 == 1,
                            l1 == 1 & l2 == 2 & l12 == 1,
                            l1 == 1 & l2 == 1 & l12 == 0,
                            l1 == 2 & l2 == 2 & l12 == 2,
                            l1 == 1 & l2 == 2 & l12 == 0,
                            l1 == 2 & l2 == 2 & l12 == 1,
                            l1 == 2 & l2 == 2 & l12 == 0))
        cat("produces", length(groups), "groups \n")
        cat("DIFF", length(groups) - c(3,2,1,4,1,2,1)[category], "\n")
        cat("ACTUAL DIFF ****", length(unique(c(Sh[r_u,s_v]))) - c(3,2,1,4,1,2,1)[category], "****\n")
      }
      
    }
  }
  
  return(Sh)
}




# fake setup (simply upload your dataset) ---------------------------------
n <- 100
d <- 16
clus <- rep(1:4, each = 4)
K <- length(unique(clus))
clus <- match(clus,1:K)

sigs <- runif(K,.1,.25)
Z <- sapply(sigs, function(s) rnorm(n,0,s))
z <- rnorm(n,0,.1)
X <- matrix(rnorm(n*d,0,.1),n,d) + Z[,clus] + z



# run algorithm -----------------------------------------------------------

n <- nrow(X)
d <- ncol(X)
p <- choose(d,2)

# compute tau and its jackknife variance in O(n*log(n))
TSh <- tautests::tau_and_jack(X)
Sh <- TSh$jack_var
Th <- TSh$tau
image(t(Th)[,d:1])
th <- Th[tautests::R2IJ(1:p)]


St <- constrainSigma(Sh, clus)
image(t(St)[,p:1])

# checkup (when all clusters are greater than 3)
length(unique(c(St))) == (K^4 + 6*K^3 + 11*K^2 + 6*K)/8



###########################################################################
# true application --------------------------------------------------------
###########################################################################

X <- as.matrix(fread("~/Downloads/data_pe.csv"))
n <- nrow(X)
d <- ncol(X)
p <- choose(d,2)

# compute tau and its jackknife variance in O(n*log(n))
TSh <- tautests::tau_and_jack(X)
Sh <- TSh$jack_var
Th <- TSh$tau
image(t(Th)[,d:1])
th <- Th[tautests::R2IJ(1:p)]

hc <- hclust(as.dist(1-abs(Th)), "ward.D2")
plot(hc)
K <- 5
clus <- cutree(hc, K)

St <- constrainSigma(Sh, clus)
image(t(St)[,p:1])

eig <- eigen(St)
plot(eig$values)
Sti <- eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
Sti2 <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)

# checkup (when all clusters are greater than 3)
length(unique(c(St))) == (K^4 + 6*K^3 + 11*K^2 + 6*K)/8



ij <- tautests::R2IJ(1:p)
A <- matrix(0,d,K)
A[cbind(1:d,clus)] <- 1  
A <- t(apply(ij, 1, function(x) colSums(A[x,])))
A_u <- unique(A)
B <- matrix(0, nrow(A), nrow(A_u))
B[cbind(1:p, apply(A, 1, function(a) which.max(colSums(a == t(A_u)))))] <- 1
B <- B[, colSums(B) != 0]
L <- ncol(B)

G <- B %*% MASS::ginv(B)
# IG <- diag(p) - G
tt <- G %*% th
Tt <- diag(d); Tt[ij] <- Tt[ij[,2:1]] <- tt
par(mfrow=c(1,2), mar=c(2,2,1,1))
image(t(Th)[,d:1], zlim = range(c(Th,Tt)))
image(t(Tt)[,d:1], zlim = range(c(Th,Tt)))

loss_E <- mahalanobis(th, tt, Sti, T)
pchisq(loss_E, p - L, lower.tail = F)

loss_M <- max(abs(Sti2 %*% (th- tt)))
mc_loss_M <- apply(mvtnorm::rmvnorm(1000, sigma = St),1,function(x) max(abs(x)))
mean(loss_M <= mc_loss_M)




