Tau.hat <- cor.fk(X)
hc <- hclust(dist(Tau.hat), method="ward.D2")
plot(hc)

oo <- hc$order
image(t(cor.fk(X)[rev(oo),oo]), zlim = c(-1,1), col = colorRampPalette(c("red","aliceblue","green"))(100))

clus <- cutree(hc,30)

Delta0 <- sapply(unique(clus), function(k){
  as.numeric(clus==k)
})

Delta <- tcrossprod(Delta0)

Tau.tilde <- (Delta %*% (Tau.hat - diag(d)) %*% Delta)/(Delta %*% (1 - diag(d)) %*% Delta)
Tt <- (t(Delta0) %*% (Tau.hat - diag(d)) %*% Delta0)/(t(Delta0) %*% (1 - diag(d)) %*% Delta0)
image(t(Tau.tilde[rev(oo),oo]), zlim = c(-1,1), col = colorRampPalette(c("red","aliceblue","green"))(100))



hc <- hclust(dist(Tau.tilde), method="ward.D2")
plot(hc)
