library(data.table)
#library(parallel)

X <- as.matrix(fread("X_NASDAQ100"))
NASDAQ100 <- fread("NASDAQ100")

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

#sapply(list.files(pattern="[.]R$", path = "C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/new code", full.names=TRUE), source)
#Theta.hat <- computeThPara(X)


#######
# running the algorithm (may take quite a while with d=106 variables)
#######

#path <- pathBuilder0Para(cor.fk(X), Theta.hat, nrow(X))
#saveRDS(path, "path_NASDAQ100_Para")

#########################################################################


path <- readRDS("path_NASDAQ100_Para")


elements <- 1:106
#a <- stoppingCriterion0Para(path$Th, path$D, path$S, elements)
#saveRDS(a, "a_NASDAQ100_Para")

a <- saveRDS("a_NASDAQ100_Para")

alpha <- .99

final <- elements[which(a < alpha)[1] - 1] 
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


par(mfrow = c(1,3), mar = c(1,1,1,1))
colfunc <- colorRampPalette(c("darkred", "darkorange", "palegoldenrod", "forestgreen", "darkblue"))
image(t(path$Th[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100), main = "The matrix Tau.hat (order = symbols)")
image(t(path$Th[new,new][d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100), main = paste0("The matrix Tau.hat (block form of ",final,")"))

D2.final <- tcrossprod(D.final,D.final)[new,new]
Tt.final <- (D2.final %*% (path$Th[new,new] - diag(d)) %*% D2.final) / (D2.final %*% (1 - diag(d)) %*% D2.final)
diag(Tt.final) <- 1

image(t(Tt.final[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(100), main = paste0("The matrix Tau.tilde ",final))

# using the data.table created in the other file, we create a very basic output of the clusters.
for(i in 1:length(new.list)){
  print(paste0(""))
  print(paste0(""))
  print(paste0(""))
  print(paste0(""))
  print(paste0("#################################"))
  print(paste0("######### cluster ", i, " ##########"))
  print(paste0("#################################"))
  print(subset(NASDAQ100[new.list[[i]], ], select = c(Symbol, Name, Sector, industry)))
  print(paste0(""))
  print(paste0("##################################################################"))
}  