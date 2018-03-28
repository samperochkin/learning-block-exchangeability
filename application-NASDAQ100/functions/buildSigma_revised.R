buildSigma <- function(Theta, tau, n){
  tau.mod <- 2*(2*n - 3) / (n*(n - 1)) * (tau + 1)^2
  Theta - tau.mod
}
