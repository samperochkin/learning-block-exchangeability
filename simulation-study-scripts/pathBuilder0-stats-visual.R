#load(paste("C:/Users/Samuel/Dropbox/HD-sim/minMSE", sep = ""))
#load(paste("C:/Users/Samuel/Dropbox/HD-sim/ThMSE", sep = ""))

# I copied to two items in the dropbox
# You should set directory to source file location (in the Session tab)
#setwd("C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/old code/Simulation_study")

load(paste("minMSE", sep = ""))
load(paste("ThMSE", sep = ""))

library(xtable)

mean(minMSE / ThMSE)

par(mfrow = c(1,1), mar = c(2,2,1,1))
hist(minMSE / ThMSE)


# These two vectors are of length 8000 = 500 x 16
# We can easily separate the cases according to n
# the scenario.grid is built such that

#     n       Tau.num
# 1   25       1
# 2   50       1
# 3   75       1
# 4  100       1
# 5   25       2
# 6   50       2
# 7   75       2


# Sooooo
# For the following matrices and figures: Tau on the x-axis, n = 25,50,75,100 on the y-axis

m <- matrix(sapply(1:16, function(i){
  1 - mean(minMSE[(i-1) * 500 + 1:500] / ThMSE[(i-1) * 500 + 1:500])
}),4,4, byrow = FALSE)



#xtable(cbind(m,m1[,1],m2[,1]))
#xtable(cbind(m1[,2],m2[,2]))



par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:16, function(i){hist(minMSE[(i-1) * 500 + 1:500] / ThMSE[(i-1) * 500 + 1:500], xlim = c(0,1), breaks = seq(0,1,.1), main = "")})

par(mfrow = c(1,1))
plot(ThMSE)
plot(minMSE)

par(mfrow = c(1,1))
plot(ThMSE, ylim = c(0,.1))
points(minMSE, col = "green")

# separated by cases

par(mfrow = c(4,4), mar = c(2,2,1,1))
sapply(1:16, function(i){
  plot(ThMSE[(i-1) * 500 + 1:500], ylim = c(0,.1))
  points(minMSE[(i-1) * 500 + 1:500], col = "green")
})


# This is weird... bu not that weird... I don't know
length(unique(ThMSE))
length(unique(minMSE))
