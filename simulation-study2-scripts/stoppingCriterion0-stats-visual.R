#load(paste("C:/Users/Samuel/Dropbox/HD-sim/minMSE", sep = ""))
#load(paste("C:/Users/Samuel/Dropbox/HD-sim/ThMSE", sep = ""))

# I copied to two items in the dropbox
# You should set directory to source file location (in the Session tab)
# setwd("C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/old code/Simulation_study")

load(paste("selAll", sep = ""))
load(paste("allMSE", sep = ""))


plot(rowMeans(selAll))
plot(rep(1:20,ncol(selAll)),c(selAll))

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:16, function(i){
  plot(rep(1:20,ncol(selAll[,(i-1) * 500 + 1:500])),c(selAll[,(i-1) * 500 + 1:500]))
})

alpha <- .05
alpha <- .25
alpha <- .5
alpha <- .75
alpha <- .95
final <- sapply(1:8000, function(i){
  rev(which(selAll[,i] > alpha))[1]
})
par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:16, function(i){
  hist(final[(i-1) * 500 + 1:500], main = "", breaks = 1:20)
})


MSEfinal <- allMSE[cbind(final,1:8000)]

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:16, function(i){
  hist(1 - (MSEfinal[(i-1) * 500 + 1:500] / allMSE[1,(i-1) * 500 + 1:500]), main = "", breaks = seq(-1,1,.1))
})

library(xtable)

MSEsfinal <- sapply(1:16, function(i){
  1 - mean(MSEfinal[(i-1) * 500 + 1:500] / allMSE[1,(i-1) * 500 + 1:500])
})


m <- matrix(MSEsfinal,4,4)

#m11 <- cbind(m,m1[,1],m2[,1])
#m22 <- cbind(m,m1[,1],m2[,1])
#m33 <- cbind(m,m1[,1],m2[,1])

# rbind the results for the other matrices before using xtable.
# xtable(cbind(m,m1[,1],m2[,1]))
# xtable(cbind(m11,m22,m33))


