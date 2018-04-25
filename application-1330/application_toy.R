library(data.table)
library(dendextend)
library(Matrix)
library(matrixcalc)
library(pcaPP)

par(mfrow=c(1,1), mar = c(2,0,0,0))
colfunc <- colorRampPalette(c("coral","aliceblue","lightgreen"))
plot(y = rep(0,1001), x = seq(-1,1,.002), col=colfunc(1001), pch=19, cex=4, xlim = c(-1, 1), ylim = c (-.005, .005), xlab = "", ylab = "", frame.plot = FALSE, yaxt = "n")


# load X -- Here I load a matrix X for a toy example, along with the true Tau to compute the squared
# error
Tau <- readRDS("application-1330/Tau_toy")
X <- readRDS("application-1330/X_toy")

n <- nrow(X)
d <- ncol(X)
p <- d * (d - 1) / 2

# compute Tau hat
Tau.hat <- cor.fk(X)

l.ij.mat <- t(combn(1:d,2))
ij.l.mat <- matrix(nrow=d,ncol=d)
ij.l.mat[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(1:p,1:p)


par(mfrow= c(1,2), mar = c(1,1,1,1))
image(t(Tau[d:1,]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))
image(t(Tau.hat[d:1,]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))


# Load functions (most of them are useless here)
sapply(list.files(pattern="[.]R$", path = "application-1330/functions", full.names=TRUE), source)

# Compute the variance of each tau -- This is later used to find a reasonable number of clusters
Theta.hat <- computeThPara(X)
Sigma.hat <- buildSigma(Theta.hat, Tau.hat[l.ij.mat], n)

# Here we do a basic hierarchical clustering. Note the distance matrix is dist(Tau.hat)
# as opposed to as.dist(1-Tau.hat^2). It makes more sense to recover blocks.
par(mfrow= c(1,1), mar = c(2,2,1,1))
hc <- hclust(dist(Tau.hat), method = "ward.D")
plot(hc)

# let us evaluate the cut producing a small number of clusters
Ks <- 20:1

a <- sapply(Ks, function(K){
  clus <- cutree(hc,K)
  cluss <- lapply(1:K, function(k){
    which(clus==k)
  })

  # Comment that for your example. That provides the true squared error.
  Tt <- matrix(sapply(1:K, function(k){
    sapply(1:K, function(kk){
      if(k == kk){
        ll <- length(cluss[[k]])
        sum(Tau.hat[cluss[[k]],cluss[[k]]] - diag(ll))/(ll*(ll-1))
      }else{
        mean(Tau.hat[cluss[[k]],cluss[[kk]]])
      }
    })
  }),K,K)
  tau.tilde <- Tt[cbind(clus[l.ij.mat[,1]],clus[l.ij.mat[,2]])]
  Tau.tilde <- diag(d)
  Tau.tilde[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(tau.tilde,tau.tilde)
  print(paste("squared error:",sum((Tau.tilde-Tau)^2)/2))
  # comment up to here
  
# Here, you should control w so that the plot (below) is monotone decreasing (or almost)...
criterion(clus, w = .75, Tau.hat, Sigma.hat, l.ij.mat, ij.l.mat)
})

a
par(mar = c(2,2,1,1))
plot(Ks,a, type="l",xlim = c(max(Ks),min(Ks)))
points(Ks,a,pch=19,cex=.75)



K.final <- 6
clus <- cutree(hc,K.final)
cluss <- lapply(1:K.final, function(k){
  which(clus==k)
})
Tt <- matrix(sapply(1:K.final, function(k){
  sapply(1:K.final, function(kk){
    if(k == kk){
      ll <- length(cluss[[k]])
      sum(Tau.hat[cluss[[k]],cluss[[k]]] - diag(ll))/(ll*(ll-1))
    }else{
      mean(Tau.hat[cluss[[k]],cluss[[kk]]])
    }
  })
}),K.final,K.final)
tau.tilde <- Tt[cbind(clus[l.ij.mat[,1]],clus[l.ij.mat[,2]])]
Tau.tilde <- diag(d)
Tau.tilde[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(tau.tilde,tau.tilde)

par(mfrow= c(1,3), mar = c(1,1,1,1))
image(t(Tau.hat[d:1,]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))
image(t(Tau.tilde[d:1,]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))
image(t(Tau[d:1,]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))

# To get a cleaner visu for your application, use
# oo <- hc$order
# par(mfrow= c(1,2), mar = c(1,1,1,1))
# image(t(Tau.hat[rev(oo),oo]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))
# image(t(Tau.tilde[rev(oo),oo]), zlim = c(-1,1), axes=FALSE, col = colfunc(100))
