# hidden uses... not ideal.
library(pcaPP)
library(MASS)
library(matrixcalc)
library(Matrix)

# source functions --------------------------------------------------------
sapply(list.files(pattern="[.]R$", path = "functions", full.names=TRUE), source)
source("functions/pathBuilder_revised.R")
source("functions/stoppingCriterion_revised.R")
source("functions/buildSigma_revised.R")
source("functions/computeTh_revised.R")
source("functions/computeTt_revised.R")



# setup -------------------------------------------------------------------
n <- 100
d <- 6
clus <- c(1,2,2,3,3,3)
K <- 3

Z1 <- rnorm(n,0,.25)
Z2 <- rnorm(n,0,.15)
Z3 <- rnorm(n,0,.25)
Z <- matrix(rnorm(n*d,0,.1),n,d)

X <- Z + cbind(Z1+Z2,Z1+Z2+Z3,Z2+Z3)[,clus]
n <- nrow(X)
d <- ncol(X)

Theta.hat <- computeTh(X)
Tau.hat <- cor.fk(X)

image(t(Tau.hat[d:1,]), zlim=c(0,1))


# run algorithm -----------------------------------------------------------
path <- pathBuilder(Tau.hat, Theta.hat, nrow(X), 1)
a <- stoppingCriterion(Tau.hat, path$D, path$S)

plot(0:(d-1),a,pch=19,xlab="d-K",ylab="critere")
lines(0:(d-1),a)
abline(h=.05,col="red")

alpha <- .05
# we select K=
which(rev(a) > alpha)[1]
# or
d-which(c(a,0) < alpha)[1]+2


# Simulation experiment ---------------------------------------------------

# wrap the previous procedure in a function called forRep and source it
source("forRep.R")

alpha <- .05
R <- 200

al <- replicate(R, {
  forRep(n)
})

dim(al) # d x R, because we have R paths of length d




# Results -----------------------------------------------------------------

# Select K
al2 <- apply(al,2,cummin)
Ks <- apply(al2,2,function(a) which(rev(a) > alpha)[1] )

# We almost always choose K=3
table(Ks)

# Reality check. Would be nice to be around alpha
mean(al[d-2,] < alpha)
alpha




# Vizu --------------------------------------------------------------------

# If hit the structure most times

xx <- (1:R)/(R+1)

plot(xx,sort(al[1,]), ylim = c(0,1.2), type = "l",col=2)
points(1/(d+1),1.05,col=2,pch=19)
text(1/(d+1),1.15,labels=paste0("K=",d),col=2,pch=19)

for(i in 2:d){
  lines(xx,sort(al[i,]), col = i+1)
  points(i/(d+1),1.05,col=i+1,pch=19)
  text(i/(d+1),1.15,labels=paste0("K=",d-i+1),col=i+1,pch=19)
}

lines(c(0,1),c(0,1), lty = 2)
