#load(paste("C:/Users/Samuel/Dropbox/HD-sim/minMSE", sep = ""))
#load(paste("C:/Users/Samuel/Dropbox/HD-sim/ThMSE", sep = ""))

# I copied to two items in the dropbox
# You should set directory to source file location (in the Session tab)
setwd("C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/old code/Simulation_study")

load(paste("allMSEHD1", sep = ""))
allMSEHD1 <- allMSEHD
load(paste("allMSEHD2", sep = ""))
allMSEHD2 <- allMSEHD
rm(allMSEHD)

minMSEHD1 <- sapply(1:4000, function(i){
  min(allMSEHD1[,i])
})
ThMSEHD1 <- allMSEHD1[1,]

minMSEHD2 <- sapply(1:4000, function(i){
  min(allMSEHD2[,i])
})
ThMSEHD2 <- allMSEHD2[1,]

mean(minMSEHD1 / ThMSEHD1)
mean(minMSEHD2 / ThMSEHD2)

hist(minMSEHD1 / ThMSEHD1)


m1 <- matrix(sapply(1:8, function(i){
  1 - mean(minMSEHD1[(i-1) * 500 + 1:500] / ThMSEHD1[(i-1) * 500 + 1:500])
}),4,2, byrow = FALSE)
m2 <- matrix(sapply(1:8, function(i){
  1 - mean(minMSEHD2[(i-1) * 500 + 1:500] / ThMSEHD2[(i-1) * 500 + 1:500])
}),4,2, byrow = FALSE)

cbind(m1,m2)

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:8, function(i){hist(minMSEHD1[(i-1) * 500 + 1:500] / ThMSEHD1[(i-1) * 500 + 1:500], main = "")})
sapply(1:8, function(i){hist(minMSEHD2[(i-1) * 500 + 1:500] / ThMSEHD2[(i-1) * 500 + 1:500], main = i + 8)})

par(mfrow = c(1,1))

plot(minMSEHD1)
plot(minMSEHD2)

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:8, function(i){
  v1 <- minMSEHD1[(i-1) * 500 + 1:500]
  v2 <- ThMSEHD1[(i-1) * 500 + 1:500]
  plot(v1, main = i, col = "green", ylim = c(min(c(v1,v2)),max(c(v1,v2))))
  points(v2)
})
sapply(1:8, function(i){
  v1 <- minMSEHD2[(i-1) * 500 + 1:500]
  v2 <- ThMSEHD2[(i-1) * 500 + 1:500]
  plot(v1, main = i, col = "green", ylim = c(min(c(v1,v2)),max(c(v1,v2))))
  points(v2)
})


