library(data.table)
library(dendextend)
#library(parallel)
library(Matrix)
library(matrixcalc)
library(pcaPP)


# load X


n <- nrow(X)
d <- ncol(X)
p <- d * (d - 1) / 2

# compute Tau hat
Tau.hat <- cor.fk(X)

# Load functions
sapply(list.files(pattern="[.]R$", path = "application-1330/functions", full.names=TRUE), source)

# That approximates the variance of each tau(X_i,X_j)
Theta.hat <- computeThPara(X)
# Since n=29 is very small, I suggest you apply shrinkage the variances (this is what the param w does)
# cutoff is the maximum number of pairs we try to merge (otherwise we have too many comparisons to try)
path <- pathBuilder(Tau.hat, Theta.hat, nrow(X), buffer = 20, cutoff = 20, w=.75)


l.ij.mat <- t(combn(1:d,2))
ij.l.mat <- matrix(nrow=d,ncol=d)
ij.l.mat[rbind(l.ij.mat,l.ij.mat[,2:1])] <- c(1:p,1:p)

ind <- which(Tau.hat[l.ij.mat] > .1)

D <- matrix(nrow=d,ncol=d)

#Sigma.hat <- buildSigma()
S <- .5*Sigma.hat + .5*mean(Sigma.hat)

D[l.ij.mat[ind,]] <- sapply(ind,function(r){
  i <- l.ij.mat[r,1]
  j <- l.ij.mat[r,2]
  sum(2*(Tau.hat[-c(i,j),i]-Tau.hat[-c(i,j),j])^2/(S[ij.l.mat[-c(i,j),i]] + S[ij.l.mat[-c(i,j),j]]))
})

D[l.ij.mat[ind,2:1]] <- D[l.ij.mat[ind,]]

D[is.na(D)] <- max(D,na.rm=TRUE)+.01

hc <- hclust(as.dist(D))




elements <- c(950:999)
# Here increase w until the selection criterion is near-monotone.
a <- stoppingCriterion(path$Th, path$clus, path$S, elements, w=1)
a
par(mar = c(2,2,1,1))
plot(elements,a, type="l")#,xlim = c(425,500))
points(elements,a)


