# This is just an example: consider a case where n = 29 and d = 13,000
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
values <- c(.2,.35,.65)
# size of clusters
d.vec <- c(90,30,1000)
d.e <- rbind(values,cumsum(c(1,d.vec[-length(d.vec)])),cumsum(d.vec))
d <- sum(d.vec)

par(mar = c(0,0,0,0))

Sigma <- constructMatrix(values,d.vec)
#image(t(Sigma[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))

Sigma <- addBlock(Sigma,d.e,1,2,.1,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,1,3,.05,needPD=TRUE)
Sigma <- addBlock(Sigma,d.e,2,3,.2,needPD=TRUE)
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
source("application-NASDAQ100/functions/pathBuilder0Para_revised.R")

Theta.hat <- computeThPara(X)


#######
# running the algorithm (may take quite a while with d=13000 variables)
#######

# Here we get an error, probably because we do not take advantage the Theta being a vector.
path <- pathBuilder0Para(cor.fk(X), Theta.hat, nrow(X), cutoff = 50)


elements <- 1:d
a <- stoppingCriterion0Para(path$Th, path$D, path$S, elements)


alpha <- .95
final <- which(a < alpha)[1]

# So the final i=
d - final + 1

#final <- 80
D.final <- path$D[[final]]


# sizes of the K = 16 clusters (notice they are all > 1)
colSums(D.final)

# So the number of distinct parameters is
K <- d - final + 1
K * (K+1) / 2

# This is a small number compared to p
p

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

#hc <- hclust(dist(Tt.small))
hc <- hclust(as.dist(1 - abs(Tt.small)))


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
  print(paste0("##################################################################"))
  if(length(new.list[[i]]) == 1){
    print(paste0("######### cluster ", i, " ################### tau.tilde = 1.000 #######"))
  }else{
    print(paste0("######### cluster ", i, " ################### tau.tilde = ",round(Tt.final[conversion[new.list[[i]][1]],conversion[new.list[[i]][2]]], digits = 3) ," #######"))
  }
  print(paste0("##################################################################"))
  #print(subset(NASDAQ100[new.list[[i]], ], select = c(Symbol, Name, Sector)))
  print(subset(NASDAQ100[new.list[[i]], ], select = c(Symbol, Name, Sector, industry)))
}  
