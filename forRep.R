forRep <- function(n){
  # setup -------------------------------------------------------------------
  d <- 6
  clus <- c(1,2,2,3,3,3)
  K <- 3
  
  Z1 <- rnorm(n,0,.25)
  Z2 <- rnorm(n,0,.15)
  Z3 <- rnorm(n,0,.25)
  Z <- matrix(rnorm(n*d,0,.1),n,d)
  
  X <- Z + cbind(Z1+Z2,Z1+Z2+Z3,Z2+Z3)[,clus]
  n <- nrow(X)
  d <- ncol(X)
  
  Theta.hat <- computeTh(X)
  Tau.hat <- cor.fk(X)
  
  # run algorithm -----------------------------------------------------------
  path <- pathBuilder(Tau.hat, Theta.hat, nrow(X), 1)
  a <- stoppingCriterion(Tau.hat, path$D, path$S)
  a
}
