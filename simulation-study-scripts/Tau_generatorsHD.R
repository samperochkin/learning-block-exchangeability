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
  #if(is.positive.definite(Sig)==FALSE & needPD){
  #  print("not positive-definite")
  #  stop(call=FALSE)
  #}else{
  #  image(t(Sig[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))
    return(Sig)
  #}
}


generator.list <- list()
true.delta <- list()

# structured HD
values <- c(.15,.3,.45,.6,.35,.8,.9,.7,.55,.2)
d.vec <- c(3,5,3,4,6,3,2,4,10,10)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

base.vector <- seq(from = .2,to = -.1, by = -.03)

set.seed(45)
inter.values <- sapply((length(values) - 1):1, function(i){
  base.vector[1:i] + rnorm(i,0,.03)
})

for(i in 1:length(inter.values)){
  for(j in 1:length(inter.values[[i]])){
    Tau <- addBlock(Tau,d.e,i,i + j,inter.values[[i]][j])
  }
}

is.positive.definite(Tau)
is.positive.definite(sin(pi / 2 * Tau))
image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

generator.list[[1]] <- Tau
true.delta[[1]] <- constructMatrix(rep(1,length(values)),d.vec)



# hard even
values <- rep(1,50)
d.vec <- rep(1,50)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

base.vector <- seq(from = .2,to = -.1, by = -.006)

set.seed(46)
inter.values <- sapply((length(values) - 1):1, function(i){
  base.vector[1:i] + rnorm(i,0,.03)
})

for(i in 1:length(inter.values)){
  for(j in 1:length(inter.values[[i]])){
    Tau <- addBlock(Tau,d.e,i,i + j,inter.values[[i]][j])
  }
}

is.positive.definite(Tau)
is.positive.definite(sin(pi / 2 * Tau))
image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))


generator.list[[2]] <- Tau
true.delta[[2]] <- constructMatrix(rep(1,length(values)),d.vec)



# case 3

values <- c(.15,.3,.45,.6,.35,.8,.9,.7,.55,.2,.15,.3,.45,.6,.35,.8,.9,.7,.55)
d.vec <- c(3,5,3,4,6,3,2,4,10,10,3,5,3,4,6,3,2,4,20)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

base.vector <- seq(from = .2,to = -.1, by = -.015)

set.seed(45)
inter.values <- sapply((length(values) - 1):1, function(i){
  base.vector[1:i] + rnorm(i,0,.03)
})

for(i in 1:length(inter.values)){
  for(j in 1:length(inter.values[[i]])){
    Tau <- addBlock(Tau,d.e,i,i + j,inter.values[[i]][j])
  }
}

is.positive.definite(Tau)
is.positive.definite(sin(pi / 2 * Tau))
image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

generator.list[[3]] <- Tau
true.delta[[3]] <- constructMatrix(rep(1,length(values)),d.vec)




# hard even
values <- rep(1,100)
d.vec <- rep(1,100)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

Tau <- constructMatrix(values,d.vec)
is.positive.definite(Tau)
image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

base.vector <- seq(from = .4,to = -.1, by = -.005)

set.seed(48)
inter.values <- sapply((length(values) - 1):1, function(i){
  base.vector[1:i] + rnorm(i,0,.1)^2
})

for(i in 1:length(inter.values)){
  for(j in 1:length(inter.values[[i]])){
    Tau <- addBlock(Tau,d.e,i,i + j,inter.values[[i]][j])
  }
}

is.positive.definite(Tau)
is.positive.definite(sin(pi / 2 * Tau))
image(t(sin(pi / 2 * Tau)[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))


generator.list[[4]] <- Tau
true.delta[[4]] <- constructMatrix(rep(1,length(values)),d.vec)

# need to remove d because not constant
rm("d")