# This is just an example: consider a case where n = 29 and d = 1330
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
values <- c(.2,.35,.65,.7,.3)
# size of clusters
d.vec <- c(200,200,200,200,200)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

par(mar = c(0,0,0,0))

Sigma <- constructMatrix(values,d.vec)
#image(t(Sigma[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

Sigma <- addBlock(Sigma,d.e,1,2,.1,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,3,.05,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,4,-.05,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,3,.2,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,4,0,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,3,4,.2,needPD=TRUE)
# Note that I should not create the whole matrix, but instead only the necessary columns (there are just 3)
# I will improve that shortly

# Bijective relation between Tau and Sigma
Tau <- 2/pi*asin(Sigma)

# total number of parameters
p <- d * (d - 1) / 2

# the data generated
n <- 29
X <- rmvnorm(n,rep(0,d),Sigma)


########
# loading the necessary functions and computing the parameter Theta
########

sapply(list.files(pattern="[.]R$", path = "simulation-study2-scripts/functions", full.names=TRUE), source)
sapply(list.files(pattern="[.]R$", path = "application-NASDAQ100/functions", full.names=TRUE), source)
source("application-1330/functions/pathBuilder.R")
source("application-1330/functions/stoppingCriterionPara.R")

Theta.hat <- computeThPara(X)


#######
# running the algorithm (may take quite a while with d=13000 variables)
#######

# Here we get an error, probably because we do not take advantage the Theta being a vector.
t1 <- Sys.time()
path <- pathBuilder(cor.fk(X), Theta.hat, nrow(X), buffer = 10, cutoff = 10, w=.6)
difftime(Sys.time(),t1)
# D2 <- tcrossprod(path$D[[485]])
# image(t(D2[d:1,]))
# Tt <- (D2 %*% (path$Th - diag(d)) %*% D2)/(D2 %*% (1 - diag(d)) %*% D2)
# diag(Tt) <- 1
# sum((path$Th-Tau)^2)/2
# sum((Tt-Tau)^2)/2

# image(t(path$Th[d:1,]), zlim = c(-1,1))
# image(t(Tt[d:1,]), zlim = c(-1,1))

#elements <- c(400,410,435:480,490)
elements <- c(950:999)
# Here increase w until the selection criterion is near-monotone.
#a <- stoppingCriterionPara(path$Th, path$D, path$S, elements, w=1)
a <- stoppingCriterion(path$Th, path$clus, path$S, elements, w=1)
a
par(mar = c(2,2,1,1))
plot(elements,a, type="l")#,xlim = c(425,500))
points(elements,a)

alpha <- .95
final <- rev(elements[which(a > alpha)])[1]


tau <- Tau[l.ij.mat]
elements <- seq(975,1000,2)
mse <- numeric(length(elements))
for(i in seq_along(elements)){
  clus <- clus.list[[i]]
  K <- max(clus)
  cluss <- lapply(1:K, function(k){
    which(clus==k)
  })
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
  print(i)
  #D <- path$D[[elements[i]]]
  #D2 <- tcrossprod(D)
  #Tt <- (D2 %*% (path$Th - diag(d)) %*% D2)/(D2 %*% (1 - diag(d)) %*% D2)
  #diag(Tt) <- 1
  tau.tilde <- Tt[cbind(clus[l.ij.mat[,1]],clus[l.ij.mat[,2]])]
  mse[i] <- sum((tau.tilde-tau)^2)
}

par(mar = c(2,2,1,1))
plot(elements,mse, type = "l")


final <- 995
clus <- clus.list[[final]]
K <- max(clus)
cluss <- lapply(1:K, function(k){
  which(clus==k)
})
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
Tau.tilde <- diag(d)
Tau.tilde[l.ij.mat] <- Tt[cbind(clus[l.ij.mat[,1]],clus[l.ij.mat[,2]])]
Tau.tilde[l.ij.mat[,2:1]] <- Tau.tilde[l.ij.mat]

# So the final i=
#D2 <- tcrossprod(D.final)
#Tt <- (D2 %*% (path$Th - diag(d)) %*% D2)/(D2 %*% (1 - diag(d)) %*% D2)
#diag(Tt) <- 1

par(mfrow= c(1,2), mar = c(1,1,1,1))
#image(t(path$Th[d:1,]), zlim = c(-1,1), xaxt = "n", yaxt = "n")
#image(t(Tau.tilde[d:1,]), zlim = c(-1,1), xaxt = "n", yaxt = "n")

sum((path$Th-Tau)^2)/2
sum((Tau.tilde-Tau)^2)/2

