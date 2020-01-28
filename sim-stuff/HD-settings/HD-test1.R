# This is just an example: consider a case where n = 70 and d = 1330
library(data.table)
library(parallel)
library(Matrix)
library(matrixcalc)
library(mvtnorm)
library(pcaPP)

#### This part you can just run ####
# They are functions I used a while ago to generate a desired matrix

# The colorpalette
colfunc <- colorRampPalette(c("coral","aliceblue","lightgreen"))
par(mar = c(0,0,0,0))
plot(y = rep(0,1001), x = seq(-1,1,.002), col=colfunc(1001), pch=19, cex=4, xlim = c(-1, 1), ylim = c (-.005, .005), xlab = "", ylab = "", frame.plot = FALSE, yaxt = "n")

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

addBlock <- function(Sig,d.e,i1,i2,par,needPD=TRUE,im=FALSE){
  Sig[d.e[2,i1]:d.e[3,i1],d.e[2,i2]:d.e[3,i2]] <- par
  Sig[d.e[2,i2]:d.e[3,i2],d.e[2,i1]:d.e[3,i1]] <- par
  if(is.positive.definite(Sig)==FALSE & needPD){
    print("not positive-definite")
    stop(call=FALSE)
  }else{
    if(im){
      image(t(Sig[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
    }
    return(Sig)
  }
}
##################################################################

# Within-cluster sigmas
values <- c(.2,.35,.65,.7,.3,.45)
# size of clusters
d.vec <- c(200,200,200,200,200,330)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

par(mar = c(0,0,0,0))

Sigma <- constructMatrix(values,d.vec)
#image(t(Sigma[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

Sigma <- addBlock(Sigma,d.e,1,2,.1,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,3,.05,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,4,-.05,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,5,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,6,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,3,.2,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,4,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,5,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,6,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,3,4,.2,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,3,5,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,3,6,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,4,5,.1,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,4,6,.1,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,5,6,.2,needPD=TRUE)
# Note that I should not create the whole matrix, but instead only the necessary columns (there are just 3)
# I will improve that shortly

# Bijective relation between Tau and Sigma
Tau <- 2/pi*asin(Sigma)

# the data generated
n <- 125
X <- rmvnorm(n,rep(0,d),Sigma)


# vectorization (i,j) -> l
l.ij.mat <- t(combn(d,2))

p <- nrow(l.ij.mat)

ij.l.mat <- matrix(nrow=d,ncol=d)
ij.l.mat[l.ij.mat] <- 1:p
ij.l.mat[l.ij.mat[,2:1]] <- ij.l.mat[l.ij.mat]

Tau.hat <- cor.fk(X)
tau.hat <- Tau.hat[l.ij.mat]

source("application-1330/functions/computeDiagThPara_revised.R")
source("application-1330/functions/buildSigma_revised.R")

Sigma.hat <-  buildSigma(computeThPara(X),tau.hat,n) 


# di <- sapply(1:nrow(l.ij.mat), function(r){
#   if(r %% 5 == 0){
#     print(r)
#   }
#   SS <- Sigma.hat[ ij.l.mat[-l.ij.mat[r,1],l.ij.mat[r,1]] ]
#   (Tau.hat[,l.ij.mat[r,1]] - Tau.hat[,l.ij.mat[r,2]])^2/4 * 
# })

Di <- dist(Tau.hat)

hc <- hclust(Di, method = "ward.D2")
plot(cutree(hc, 6))

# At K=5 we do not pass the premiliminary test, but at K=6 we do.
K <- 6
oo <- cutree(hc, K)

par(mar = c(2,2,1,1))
plot(diff(hc$height))

kk.mat <- matrix(oo[c(l.ij.mat)],ncol=2)
kk.mat <- t(apply(kk.mat,1,sort))

u.kk.mat <- unique(kk.mat)
          
blocks <- sapply(1:nrow(u.kk.mat), function(r){
  which(kk.mat[,1]==u.kk.mat[r,1] & kk.mat[,2]==u.kk.mat[r,2])
})

tau.tilde <- tau.hat
for(k in 1:length(blocks)){
  tau.tilde[blocks[[k]]] <- mean(tau.hat[blocks[[k]]])
}


# Rough value of the test
al <- 1-pchisq(sum((tau.hat-tau.tilde)^2/Sigma.hat), p - K*(K-1)/2)
al
sum((tau.hat-tau.tilde)^2/Sigma.hat)


alpha <- .05

# Test as if they were all independent (using Sidak)
alpha <- .05
#al > (1-alpha)^(1/p)



b <- blocks[[19]]
plot(Sigma.hat[b])

scores <- do.call(rbind, lapply(blocks, function(b){
  cbind(b,(tau.hat[b]-tau.tilde[b])^2/mean(Sigma.hat[b]))
}))

scores <- scores[order(scores[,1]),2]

# But they aren't. So we divide them in more or less independent clusters.
# We compute tests using the full Sigma in the clusters.

# test.groups <- 1:d #sample(hc$order)
# test.groups <- split(test.groups, ceiling(seq_along(test.groups)/20))

# source("functions/computeTh.R")
# source("functions/computeTt.R")

# Sigmas <- lapply(test.groups, function(tg){
#   Th <- computeTh(X[,tg])
# 
#   Delta <- t(sapply(tg, function(i){
#     vec <- rep(0,K)
#     vec[oo[i]] <- 1
#     vec
#   }))
#   
# 
#   Tt <- computeTt()
# })

test.groups <- sample(hc$order)
# test.groups <- 1:d
test.groups <- split(test.groups, ceiling(seq_along(test.groups)/200))

blocks <- lapply(test.groups, function(tg){
 ij.l.mat[t(combn(tg,2))]
})

al <- sapply(blocks, function(b){
  #sum((tau.hat[b]-tau.tilde[b])^2/Sigma.hat[b])
  #score <- sum((tau.hat[b]-tau.tilde[b])^2/Sigma.hat[b])
  #print(length(b) - length(unique(tau.tilde[b])))
  #1-pchisq(score, length(b) - length(unique(tau.tilde[b])))
  1-pchisq(sum(scores[b]), length(b))# - length(unique(tau.tilde[b])))
})
taus <- sapply(blocks, function(b){
  mean(tau.tilde[b])
})


plot(al)
plot(sort(al))
plot(taus,al)
hist(al,xlim=c(0,1), breaks=seq(0,1,.05))
sum(al < .05)/length(al)

# That provides many tests that we all want true.
# If these tests are (truly) independent enough, we can use Sidak's correction.
prod(al < 1-(1-.05)^(1/K))

# Now here K does not has to match the number of clusters hypothesied about...
# We could go even more hardcore and consider the number of test to grow too.



b <- blocks[[10]]

cbind(tau.hat[b],tau.tilde[b],Sigma.hat[b])

cbind((tau.hat[1:100]-tau.tilde[1:100])^2,Sigma.hat[1:100])

# Once 


hist(tau.hat[b], probability=TRUE)
hist(tau.hat[b]-tau.tilde[b], probability=TRUE)
hist((tau.hat[b]-tau.tilde[b])/sqrt(mean(Sigma.hat[b])), probability=TRUE)
lines(seq(-4,4,.01),dnorm(seq(-4,4,.01),0,1))

mean(Sigma.hat[b])
mean((tau.hat[b]-tau.tilde[b])^2)


# Create bootstrap samples of the variables (pick according to importance...)
# Perform test for each subset
#


max(Tau.hat[test.groups[[1]],test.groups[[2]]])
range(Tau.hat[test.groups[[1]],test.groups[[2]]])
mean(Tau.hat[test.groups[[1]],test.groups[[2]]])
    
max(Tau(test.groups[[i]],test.groups[[j]]))
mean(Tau(test.groups[[i]],test.groups[[j]]))
