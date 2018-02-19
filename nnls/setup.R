library(data.table)
library(parallel)
library(Matrix)
library(MASS)
library(pcaPP)
library(matrixcalc)

#colfunc <- colorRampPalette(c("coral","aliceblue","lightgreen"))
#colfunc <- colorRampPalette(c("red","lightyellow","green"))
colfunc <- colorRampPalette(c("darkred","lightyellow","forestgreen"))
plot(y = rep(0,1001), x = seq(-1,1,.002), col=colfunc(1001), pch=19, cex=4, xlim = c(-1, 1), ylim = c (-.005, .005), xlab = "", ylab = "", frame.plot = FALSE, yaxt = "n")
#plot(x = rep(0,501), y = seq(1,0,-.002), col=colfunc2(501), pch=19, cex=4, ylim = c(-1.1, 1.1), xlim = c (-.005, .005), xlab = "", ylab = "", frame.plot = FALSE, xaxt = "n")
#points(x = rep(0,501), y = seq(-1,0,.002), col=colfunc3(501), pch=19, cex=4)


constructMatrix <- function(values,d.vec){
  k <- length(values)
  d <- sum(d.vec)
  first <- matrix(c(rep(values[1],d.vec[1]),rep(0,d-d.vec[1])),byrow=TRUE,ncol=d,nrow=d.vec[1])
  middle <- sapply(2:(k-1),function(j){matrix(c(rep(0,sum(d.vec[1:(j-1)])),rep(values[j],d.vec[j]),rep(0,sum(d.vec[(j+1):k]))),byrow=TRUE,ncol=d,nrow=d.vec[j])},simplify = FALSE)
  last <- matrix(c(rep(0,sum(d.vec[1:(k-1)])),rep(values[k],d.vec[k])),byrow=TRUE,ncol=d,nrow=d.vec[k])
  
  M <- rbind(first,do.call(rbind,middle),last)
  diag(M) <- 1
  M
}

addBlock <- function(Sig,d.e,i1,i2,par,needPD = TRUE){
  Sig[d.e[2,i1]:d.e[3,i1],d.e[2,i2]:d.e[3,i2]] <- par
  Sig[d.e[2,i2]:d.e[3,i2],d.e[2,i1]:d.e[3,i1]] <- par
  if(is.positive.definite(Sig)==FALSE & needPD){
    print("not positive-definite")
    stop(call=FALSE)
  }else{
    image(t(Sig[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))
    return(Sig)
  }
}


values <- c(.4,.4,.45,.5,.2,1)
d.vec <- c(4,3,4,2,5,1)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
par(mar = c(0,0,0,0))

Tau <- addBlock(Tau,d.e,1,2,.25,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,3,.15,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,3,.2,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,4,.05,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,4,.05,needPD = FALSE)
Tau <- addBlock(Tau,d.e,3,4,.05,needPD = FALSE)




values <- c(.4,.4,.45,.5,.2,1,.8,.7)
d.vec <- c(4,3,4,2,5,1,2,3)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
par(mar = c(0,0,0,0))

Tau <- addBlock(Tau,d.e,1,2,.25,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,3,.25,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,3,.15,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,4,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,4,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,3,4,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,5,-.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,5,-.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,3,5,-.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,4,5,-.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,6,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,6,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,3,6,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,4,6,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,5,6,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,7,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,7,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,3,7,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,4,7,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,5,7,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,1,8,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,2,8,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,3,8,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,4,8,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,5,8,.1,needPD = FALSE)
Tau <- addBlock(Tau,d.e,6,7,.5,needPD = FALSE)
Tau <- addBlock(Tau,d.e,6,8,.25,needPD = FALSE)
Tau <- addBlock(Tau,d.e,7,8,.35,needPD = FALSE)

Sig <- sin(pi / 2 * Tau)


n <- 80
X <- mvtnorm::rmvnorm(n, rep(0,d), sigma = Sig)

#sapply(list.files(pattern="[.]R$", path = "simulation-study2-scripts/functions", full.names=TRUE), source)
#Theta.hat <- computeTh(X)
#sapply(list.files(pattern="[.]R$", path = "application-NASDAQ100/functions", full.names=TRUE), source)
#Theta.hat <- computeThPara(X)
l.ij.mat <- t(combn(1:d,2))

sapply(list.files(pattern="[.]R$", path = "functions", full.names=TRUE), source)
Tau.hat <- cor.fk(X)
tau.hat <- Tau.hat[l.ij.mat]
Theta.hat <- computeTh(X)
Sigma.hat <- buildSigma(Theta.hat, tau.hat, n)

Delta <- t(unique(constructMatrix(rep(1,4),c(4,3,4,2))))
delta <- tcrossprod(Delta)[l.ij.mat]

par(mar = c(1,1,1,1), mfrow = c(1,3))
image(t(Delta[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))


Gamma <- tcrossprod(delta,delta)




library(nnls)

p <- d*(d-1)/2
delta.hat <- rep(0,p)


K <- 4
lambda = 1000

A <- matrix(tau.hat/diag(Sigma.hat),p,p, byrow = TRUE)
A <- rbind(A,rep(1,d))

sol <- nnnpls(A,b = c(tau.hat/diag(Sigma.hat),K*(K+1)/2*lambda) ,con = c(rep(1,p),1))
delta.hat <- sol$x + delta.hat
A = delta.hat %*% t(tau.hat)
A <- rbind(A,rep(1,d))
#



n <- 100
X <- mvtnorm::rmvnorm(n, rep(0,d), sigma = Sig)

l.ij.mat <- t(combn(1:d,2))
sapply(list.files(pattern="[.]R$", path = "functions", full.names=TRUE), source)
Tau.hat <- cor.fk(X)
tau.hat <- Tau.hat[l.ij.mat]
Theta.hat <- computeTh(X)
Sigma.hat <- buildSigma(Theta.hat, tau.hat, n)

U <- apply(X,2,rank)
U <- t(apply(U,1,rank))
Rho.hat <- cor.fk(U)
Tau.hat <- cor.fk(X)
par(mfrow = c(1,3), mar = c(0,1,0,1))
image(t(Rho.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

plot(hclust(as.dist(1-Rho.hat)))
plot(hclust(as.dist(1-Tau.hat)))

n <- 50
X <- mvtnorm::rmvnorm(n, rep(0,d), sigma = Sig)

X2 <- X
X2[,14:18] <- -X2[,14:18]
U <- apply(X2,2,rank)
U <- t(apply(U,1,rank))
Rho.hat <- cor.fk(U)
Tau.hat <- cor.fk(X2)
par(mfrow = c(1,3), mar = c(0,1,0,1))
image(t(Rho.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
plot(hclust(as.dist(1-Rho.hat)))
plot(hclust(as.dist(1-Tau.hat)))
plot(hclust(dist(t(X2))))
plot(hclust(dist(t(X))))

Rho.hat[1:13,14:d]
Rho.hat[1:11,1:11]

plot(hclust(as.dist(1-Rho.hat)))
plot(hclust(as.dist(1-Tau.hat)))



# sparse pca

library(sparsePCA)


n <- 50
X <- mvtnorm::rmvnorm(n, rep(0,d), sigma = Sig)

l.ij.mat <- t(combn(1:d,2))
sapply(list.files(pattern="[.]R$", path = "functions", full.names=TRUE), source)
Tau.hat <- cor.fk(X)
tau.hat <- Tau.hat[l.ij.mat]
Theta.hat <- computeTh(X)
Sigma.hat <- buildSigma(Theta.hat, tau.hat, n)


m <- 6
sparsePCA(Tau.hat,m,rep(.1,m))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))


p <- d*(d-1)/2

ij.l.mat <- matrix(NA,d,d)
ij.l.mat[l.ij.mat] <- 1:p
ij.l.mat[l.ij.mat[,2:1]] <- 1:p

a.vec <- sapply(1:p, function(r){
  ii <- l.ij.mat[r,1]
  jj <- l.ij.mat[r,2]
  Tau.tilde <- (Tau.hat[ii,-c(ii,jj)] + Tau.hat[jj,-c(ii,jj)])/2
  #res <- (Tau.hat[ii,-c(ii,jj)] - Tau.tilde)^2/diag(Sigma.hat[ij.l.mat[ii,-c(ii,jj)],ij.l.mat[ii,-c(ii,jj)]]) +
  #  (Tau.hat[jj,-c(ii,jj)] - Tau.tilde)^2/diag(Sigma.hat[ij.l.mat[jj,-c(ii,jj)],ij.l.mat[jj,-c(ii,jj)]])
  res <- t(Tau.hat[ii,-c(ii,jj)] - Tau.tilde) %*% solve(Sigma.hat[ij.l.mat[ii,-c(ii,jj)],ij.l.mat[ii,-c(ii,jj)]]) %*% (Tau.hat[ii,-c(ii,jj)] - Tau.tilde) +
    t(Tau.hat[jj,-c(ii,jj)] - Tau.tilde) %*% solve(Sigma.hat[ij.l.mat[jj,-c(ii,jj)],ij.l.mat[jj,-c(ii,jj)]]) %*% (Tau.hat[jj,-c(ii,jj)] - Tau.tilde)
  1 - pchisq(sum(res), d-2)
})

A <- diag(d)
A[l.ij.mat] <- a.vec
A[l.ij.mat[,2:1]] <- a.vec
A
image(t(A[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

par(mfrow = c(1,1))

m <- 8
pcad <- sparsePCA(A,m,rep(.75,m))
pcad$Y
image(t(pcad$Y[d:1,]), axes=FALSE, col=colfunc(100))

plot(hclust(dist(pcad$Y), method = "complete"))
plot(hclust(dist(pcad$Y), method = "average"))


Delta <- t(unique(constructMatrix(rep(1,length(values)),d.vec)))
delta <- tcrossprod(Delta,Delta)[l.ij.mat]


w <- .25
Sigma.tilde <- w* diag(diag(Sigma.hat)) + (1-w)*Sigma.hat
is.positive.definite(Sigma.tilde)

Si <- solve(Sigma.tilde)
A <- Si %*% tau.hat

Delta.hat <- diag(d)
tau.tilde <- tau.hat

#delta <- rep(1/p,p)
delta <- rep(0,p)
Delta.hat[l.ij.mat[which.max(A),1],l.ij.mat[which.max(A),2]] <- 1
Delta.hat[l.ij.mat[which.max(A),2],l.ij.mat[which.max(A),1]] <- 1
delta <- Delta.hat[l.ij.mat]
tau.tilde <- t((Delta.hat %*% (Tau.hat - diag(d)) %*% Delta.hat)/(Delta.hat %*% (1 - diag(d)) %*% Delta.hat))[l.ij.mat]

der <- sapply(1:p, function(r){
  #if(Delta.hat[l.ij.mat[r,1],l.ij.mat[r,2]] == 1){
  #  return(0)
  #}
  er <- rep(0,p)
  er[r] <- 1
  
  # Remember that we want to maximize. So we look for positice derivative.
  #2*t(tau.hat) %*% delta %*% t(er) %*% Si %*% tau.hat
  #A[r]*t(tau.hat) %*% delta
  A[r]*t(tau.tilde) %*% delta
})


plot(der)

#l.ij.mat[which.max(der),]
l.ij.mat[order(der, decreasing = TRUE)[1:10],]
tcrossprod(Delta)[l.ij.mat[order(der, decreasing = TRUE),]]

Delta.hat[l.ij.mat[which.max(der),1],l.ij.mat[which.max(der),2]] <- Delta.hat[l.ij.mat[which.max(der),1],l.ij.mat[which.max(der),2]] + 1
Delta.hat[l.ij.mat[which.max(der),2],l.ij.mat[which.max(der),1]] <-Delta.hat[l.ij.mat[which.max(der),2],l.ij.mat[which.max(der),1]] + 1
delta <- Delta.hat[l.ij.mat]
tau.tilde <- t((Delta.hat %*% (Tau.hat - diag(d)) %*% Delta.hat)/(Delta.hat %*% (1 - diag(d)) %*% Delta.hat))[l.ij.mat]




delta[which.max(der)] <- delta[which.max(der)] + .01
delta

der.mat <- matrix(0,d,d)
der.mat[l.ij.mat] <- der
der.mat[l.ij.mat[,2:1]] <- der

which.max(rowSums(der.mat))
image(t(der.mat[d:1,]), axes=FALSE, col=colfunc(100))



# update delta + epsilon
# the loop
