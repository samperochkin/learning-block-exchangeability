# buildSigma <- function(Theta, tau, n){
#   p <- length(tau)
#   tau.mod <- matrix(tau, p, p) + 1
#   tau.mod <- 2*(2*n - 3) / (n*(n - 1)) * tau.mod * t(tau.mod)
#   tau.mod <- (tau.mod + t(tau.mod)) / 2
#   Theta - tau.mod
# }


buildSigma <- function(Theta, tau, n, diagOnly=TRUE){
  if(diagOnly) tau.mod <- 2*(2*n - 3) / (n*(n - 1)) * (tau + 1)^2
  else tau.mod <- 2*(2*n - 3) / (n*(n - 1)) * outer(tau+1,tau+1,"*")
  
  Theta - tau.mod
}
