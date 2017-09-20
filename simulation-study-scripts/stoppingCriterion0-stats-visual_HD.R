#load(paste("C:/Users/Samuel/Dropbox/HD-sim/minMSE", sep = ""))
#load(paste("C:/Users/Samuel/Dropbox/HD-sim/ThMSE", sep = ""))

# I copied to two items in the dropbox
# You should set directory to source file location (in the Session tab)
# setwd("C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/old code/Simulation_study")

load(paste("selAllHD1", sep = ""))
selAll1 <- selAll
load(paste("selAllHD2", sep = ""))
selAll2 <- selAll
rm(selAll)

load(paste("allMSEHD1", sep = ""))
allMSE1 <- allMSEHD
load(paste("allMSEHD2", sep = ""))
allMSE2 <- allMSEHD
rm(allMSEHD)



plot(rowMeans(selAll1))
plot(rowMeans(selAll2))
plot(rep(1:50,ncol(selAll1)),c(selAll1))
plot(rep(1:100,ncol(selAll2)),c(selAll2))

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:8, function(i){
  plot(rep(1:50,ncol(selAll1[,(i-1) * 500 + 1:500])),c(selAll1[,(i-1) * 500 + 1:500]))
})
sapply(1:8, function(i){
  plot(rep(1:100,ncol(selAll2[,(i-1) * 500 + 1:500])),c(selAll2[,(i-1) * 500 + 1:500]))
})

alpha <- .05
#alpha <- .25
#alpha <- .5
#alpha <- .75
#alpha <- .95
final1 <- sapply(1:4000, function(i){
  rev(which(selAll1[,i] > alpha))[1]
})
final2 <- sapply(1:4000, function(i){
  rev(which(selAll2[,i] > alpha))[1]
})
par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:8, function(i){
  hist(final1[(i-1) * 500 + 1:500], main = "", breaks = 1:50)
})
sapply(1:8, function(i){
  hist(final2[(i-1) * 500 + 1:500], main = "", breaks = 1:100)
})


par(mfrow = c(1,1), mar = c(2,2,1,1))
sapply(1:8, function(i){
  lines(c(selAll1[i,]))
  lines(c(selAll2[i,]))
})


MSEfinal1 <- allMSE1[cbind(final1,1:4000)]
MSEfinal2 <- allMSE2[cbind(final2,1:4000)]

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:8, function(i){
  hist(MSEfinal1[(i-1) * 500 + 1:500] / allMSE1[1,(i-1) * 500 + 1:500], main = "", breaks = seq(0,2,.1))
})
sapply(1:8, function(i){
  hist(MSEfinal2[(i-1) * 500 + 1:500] / allMSE2[1,(i-1) * 500 + 1:500], main = "", breaks = seq(0,2,.1))
})



m1 <- matrix(sapply(1:8, function(i){
  1 - mean(MSEfinal1[(i-1) * 500 + 1:500] / allMSE1[1,(i-1) * 500 + 1:500])
}),4,2)
m2 <- matrix(sapply(1:8, function(i){
  1 - mean(MSEfinal2[(i-1) * 500 + 1:500] / allMSE2[1,(i-1) * 500 + 1:500])
}),4,2)

#xtable(cbind(m1[,2],m2[,2]))


