library(MASS)
library(pcaPP)
library(matrixcalc)

colfunc <- colorRampPalette(c("darkred", "darkorange", "palegoldenrod", "forestgreen", "darkblue"))
par(mfrow = c(1,1), mar = c(2.5,0,0,0))
plot(y = rep(0,1001), x = c(c(sapply(500:1, function(i){c(i,-i)})) * .002,0), col=colfunc(1001)[c(c(sapply(500:1, function(i){c(i,-i)})),0) + 501], pch=19, cex=4, xlim = c(-1, 1), ylim = c (-.005, .005), xlab = "", ylab = "", frame.plot = FALSE, yaxt = "n")

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
    image(t(Sig[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
    return(Sig)
  }
}


generator.list <- list()
true.delta <- list()

# easy even
values <- c(.15,.3,.45,.6)
d.vec <- c(5,5,5,5)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

Tau <- addBlock(Tau,d.e,1,2,.1)
Tau <- addBlock(Tau,d.e,1,3,.05)
Tau <- addBlock(Tau,d.e,1,4,0)
Tau <- addBlock(Tau,d.e,2,3,.2)
Tau <- addBlock(Tau,d.e,2,4,.1)
Tau <- addBlock(Tau,d.e,3,4,.3)

image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

generator.list[[1]] <- Tau
true.delta[[1]] <- constructMatrix(rep(1,length(values)),d.vec)

# hard even
values <- c(.1,.2,.25,.3)
d.vec <- c(5,5,5,5)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

Tau <- addBlock(Tau,d.e,1,2,.05)
Tau <- addBlock(Tau,d.e,1,3,0)
Tau <- addBlock(Tau,d.e,1,4,-.05)
Tau <- addBlock(Tau,d.e,2,3,.1)
Tau <- addBlock(Tau,d.e,2,4,.15)
Tau <- addBlock(Tau,d.e,3,4,.2)

image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

generator.list[[2]] <- Tau
true.delta[[2]] <- constructMatrix(rep(1,length(values)),d.vec)


# hard uneven
values <- c(.1,.2,.3,0,.5,.4,.25,.05)
d.vec <- c(4,2,2,1,3,3,2,3)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

Tau <- addBlock(Tau,d.e,2,3,.15)
Tau <- addBlock(Tau,d.e,2,4,.05)
Tau <- addBlock(Tau,d.e,2,5,.1)
Tau <- addBlock(Tau,d.e,2,6,.05)
Tau <- addBlock(Tau,d.e,3,4,.25)
Tau <- addBlock(Tau,d.e,3,5,.1)
Tau <- addBlock(Tau,d.e,3,6,.15)
Tau <- addBlock(Tau,d.e,4,5,.35)
Tau <- addBlock(Tau,d.e,4,6,.25)
Tau <- addBlock(Tau,d.e,5,6,.3)
Tau <- addBlock(Tau,d.e,6,7,.15)
Tau <- addBlock(Tau,d.e,6,8,-.1)
Tau <- addBlock(Tau,d.e,7,8,-.05)

image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

generator.list[[3]] <- Tau
true.delta[[3]] <- constructMatrix(rep(1,length(values)),d.vec)


# very hard uneven
values <- c(.05,.1,.15,.2,.25,.3,.35,.4)
d.vec <- c(4,3,2,1,2,3,2,3)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

Tau <- addBlock(Tau,d.e,2,3,0)
Tau <- addBlock(Tau,d.e,2,4,.01)
Tau <- addBlock(Tau,d.e,2,5,.02)
Tau <- addBlock(Tau,d.e,2,6,.03)
Tau <- addBlock(Tau,d.e,2,7,.04)
Tau <- addBlock(Tau,d.e,2,8,.05)
Tau <- addBlock(Tau,d.e,3,4,.11)
Tau <- addBlock(Tau,d.e,3,5,.12)
Tau <- addBlock(Tau,d.e,3,6,.13)
Tau <- addBlock(Tau,d.e,3,7,.14)
Tau <- addBlock(Tau,d.e,3,8,.15)
Tau <- addBlock(Tau,d.e,4,5,.17)
Tau <- addBlock(Tau,d.e,4,6,.18)
Tau <- addBlock(Tau,d.e,4,7,.19)
Tau <- addBlock(Tau,d.e,4,8,.2)
Tau <- addBlock(Tau,d.e,5,6,.23)
Tau <- addBlock(Tau,d.e,5,7,.24)
Tau <- addBlock(Tau,d.e,5,8,.25)
Tau <- addBlock(Tau,d.e,6,7,.28)
Tau <- addBlock(Tau,d.e,6,8,.29)
Tau <- addBlock(Tau,d.e,7,8,.3)

image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

generator.list[[4]] <- Tau
true.delta[[4]] <- constructMatrix(rep(1,length(values)),d.vec)


#i <- 4
#image(t(generator.list[[i]][d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
#image(t(true.delta[[i]][d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
