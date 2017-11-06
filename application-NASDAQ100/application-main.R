library(data.table)
library(parallel)
library(Matrix)

X <- as.matrix(fread("application-NASDAQ100/X_NASDAQ100"))
NASDAQ100 <- fread("application-NASDAQ100/NASDAQ100")

full.names <- as.vector(NASDAQ100$Name)
sectors <- as.vector(NASDAQ100$Sector)
industry <- as.vector(NASDAQ100$industry)

symbols <- colnames(X)
n <- nrow(X)
d <- ncol(X)

########## No need to run ######### result available below ##############

########
# loading the necessary functions and computing the parameter Theta
########

#sapply(list.files(pattern="[.]R$", path = "simulation-study2-scripts/functions", full.names=TRUE), source)
#sapply(list.files(pattern="[.]R$", path = "application-NASDAQ100/functions", full.names=TRUE), source)
#Theta.hat <- computeThPara(X)


#######
# running the algorithm (may take quite a while with d=106 variables)
#######

#path <- pathBuilder0Para(cor.fk(X), Theta.hat, nrow(X))
#saveRDS(path, "path_NASDAQ100_Para")

#########################################################################


path <- readRDS("application-NASDAQ100/path_NASDAQ100_Para")


#elements <- 1:106
#a <- stoppingCriterion0Para(path$Th, path$D, path$S, elements)
#saveRDS(a, "a_NASDAQ100_Para")

a <- readRDS("a_NASDAQ100_Para")

par(mfrow = c(1,1), mar = c(2,2,1,1))
plot(a[87:97], pch = 19, type = "l", xaxt = "n")
axis(1, at=seq(1,d-87+1,2), labels=seq(d-87+1,1,-2))

points(x = 4:6 + 1, y = a[91:93], pch = 19)
points(x = 0:10 + 1, y = a[87:97], pch = 19)


alpha <- .95
final <- which(a < alpha)[1]

# for some very intuitive groups, run
#final <- 80

D.final <- path$D[[final]]


######
# produce the list of clusters
######

new.list <- sapply(1:ncol(D.final), function(i){
  which(D.final[,i] == 1)
})


#####
# compute the values of Tau tilde
#####

pairs <- t(combn(1:length(new.list),2))
tt.small <- sapply(1:nrow(pairs), function(i){
  j1 <- pairs[i,1]
  j2 <- pairs[i,2]
  mean(path$Th[new.list[[j1]],new.list[[j2]]])
})

#####
# construct a mini version of it with reduced blocks
#####

Tt.small <- matrix(nrow = length(new.list), ncol = length(new.list))
Tt.small[rbind(pairs,pairs[,2:1])] <- tt.small


######
# fill the diagonal with zeros...
######

diag(Tt.small) <- 0

######
# simple hierarchical clustering. One could try to modify the (main) algorithm of the paper
# so that it can be reused here. This would produce an even more informative clustering, 
# more like embbedded clusterings. See my future work for some work of this
# kind. I don't exactly do as suggested here however...
# Anywas...
# Choosing the diagonal of Tt.small to apply hclust is not so obvious. The hclust algorithm should
# h c should be a bit modified to be adequate in this case. It still gives 
# something reasonable.
######

hc <- hclust(dist(Tt.small))


# and keeping the ordering produced

new <- unlist(sapply(hc$order, function(k){
  new.list[[k]]
}))

#par(mfrow = c(1,1), mar = c(0,0,0,0))
par(mfrow = c(1,3), mar = c(1,1,1,1))
colfunc <- colorRampPalette(c("darkred", "darkorange", "palegoldenrod", "forestgreen", "darkblue"))
image(t(path$Th[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100)) # original Tau.hat
image(t(path$Th[new,new][d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100)) # newly ordered Tau.hat

D2.final <- tcrossprod(D.final,D.final)[new,new]
Tt.final <- (D2.final %*% (path$Th[new,new] - diag(d)) %*% D2.final) / (D2.final %*% (1 - diag(d)) %*% D2.final)
diag(Tt.final) <- 1

image(t(Tt.final[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100))

# to pick the right entries in Tt.final
conversion <- sapply(1:d, function(i){
  which(new == i)
})

# using the data.table created in the other file, we create a very basic output of the clusters.
for(i in 1:length(new.list)){
  print(paste0(""))
  print(paste0(""))
  print(paste0(""))
  print(paste0(""))
  print(paste0("##################################################################"))
  if(length(new.list[[i]]) == 1){
    print(paste0("######### cluster ", i, " ################### tau.tilde = 1.000 #######"))
  }else{
    print(paste0("######### cluster ", i, " ################### tau.tilde = ",round(Tt.final[conversion[new.list[[i]][1]],conversion[new.list[[i]][2]]], digits = 3) ," #######"))
  }
  print(paste0("##################################################################"))
  print(subset(NASDAQ100[new.list[[i]], ], select = c(Symbol, Name, Sector, industry)))
  print(paste0(""))
  print(paste0("###################################################################################################"))
}  
