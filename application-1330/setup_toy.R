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
n <- 29
X <- rmvnorm(n,rep(0,d),Sigma)

saveRDS(Tau,"application-1330/Tau_toy")
saveRDS(X,"application-1330/X_toy")
