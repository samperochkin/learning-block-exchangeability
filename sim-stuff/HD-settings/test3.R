dim(X)

n <- 150
X <- rmvnorm(n,rep(0,d),Sigma)

Delta <- constructMatrix(rep(1,length(d.vec)),d.vec)

Th1 <- cor.fk(X[1:125,])
Th2 <- cor.fk(X[126:150,])
Tt2 <- (Delta %*% (Th2-diag(d)) %*% Delta)/(Delta %*% (1-diag(d)) %*% Delta)

sum((Tt2-Th1)^2)/sum((Th2-Th1)^2)
sum((Tt2-Tau)^2)/sum((Th2-Tau)^2)






n <- 200
X <- rmvnorm(n,rep(0,d),Sigma)

Delta <- constructMatrix(rep(1,length(d.vec)),d.vec)

Ths <- lapply(1:2,  function(k){
  print(k)
  cor.fk(X[1:100 + (k-1)*100,])
})

# Tts <- lapply(Ths,  function(Th){
#   (Delta %*% (Th-diag(d)) %*% Delta)/(Delta %*% (1-diag(d)) %*% Delta)
# })
  
# pairs <- t(combn(1:50,2))
# 
# par(mar = c(2,2,1,1))
# hist(sapply(1:nrow(pairs), function(k){
#   sum((Ths[[pairs[k,1]]]-Tts[[pairs[k,2]]])^2)/sum((Ths[[pairs[k,1]]]-Ths[[pairs[k,2]]])^2)
# }), breaks=10)


hc <- hclust(dist(Ths[[1]]),method="ward.D2")

Ds <- lapply(1:100, function(k){
  print(k)
  clus <- cutree(hc,k)
  D <- sapply(1:k, function(i){
    vec <- rep(0,d)
    vec[which(clus==i)] <- 1
    vec
  })
  tcrossprod(D)
})


Tts <- lapply(Ds,  function(D){
  print("hey")
  (D %*% (Ths[[1]]-diag(d)) %*% D)/(D %*% (1-diag(d)) %*% D)
})
Tts <- lapply(Tts, function(Tt){
  diag(Tt) <- 1
  Tt
})

# clus <- cutree(hc,6)
# D <- sapply(1:6, function(i){
#   vec <- rep(0,d)
#   vec[which(clus==i)] <- 1
#   vec
# })
# D <- tcrossprod(D)
# Tts.hc[[6]] <- (D %*% (Ths[[1]]-diag(d)) %*% D)/(D %*% (1-diag(d)) %*% D)



par(mar = c(2,2,1,1))
SS2 <- sapply(seq_along(Tts[-1]), function(i){
  sum((Ths[[2]]-Tts[[1+i]])^2)
})
SS1 <- sapply(seq_along(Tts[-1]), function(i){
  sum((Tau-Tts[[1+i]])^2)
})

plot(SS1)
which.min(SS2)

plot(SS2)
which.min(SS2)




