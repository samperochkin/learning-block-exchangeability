# data and basic quantities
n <- 100
d <- d0 <- 8
p <- choose(d,2)
clus <- rep(1:2, each = d/2)
K <- length(unique(clus))
ij <- do.call("rbind", c(
  lapply(1:K, function(k) t(combn((k-1)*d/K + 1:(d/K),2))),
  lapply(1:(K-1), function(k){
    a1 <- (k-1)*d/K; a2 <- 1:(d/K)
    as.matrix(expand.grid(a1+a2, (a1+d/K+1):d))
  })
))

# sigs <- runif(K,.1,.25)
sigs <- rep(.1,K)
Z <- sapply(sigs, function(s) rnorm(n,0,s))
z <- rnorm(n,0,.5)
X <- matrix(rnorm(n*d,0,.1),n,d) + Z[,clus] + z
TSh <- tautests::tau_and_jack(X, ijs = ij)
Sh <- TSh$jack_var
Th <- TSh$tau
th <- Th[ij]
St <- constrainSigma(Sh, clus, ij)


# 
a <- 2*choose(d/2,2)
b <- p - a
c <- mean(St[1,1:a])
d <- mean(St[(a+1),1:a])
e <- mean(St[(a+1),(a+1):p])

sq <- sqrt(a^2*c^2 - 2*a*b*c*e + 4*a*b*d^2 + b^2*e^2)
x = (-sq + a*c - 2*a*d - b*e)/(2*a*(a*c - a*d + b*d - b*e))
y = (sq + a*c - 2*a*d - b*e)/(2*b*(a*c - a*d + b*d - b*e))
z = (-sq + a*c + b*e)/2

v <- c(rep(x, a), rep(y, b))
St %*% v - z * v



eig <- eigen(St)
ev <- round(eig$values, 5)
ev_u <- unique(ev)
ev_u

1/2 * (sqrt(a^2*c^2 - 2*a* b* c* e + 4* a* b* d^2 + b^2* e^2) + a *c + b* e)
1/2 *(-sqrt(a^2*c^2 - 2*a* b* c* e + 4* a* b* d^2 + b^2* e^2) + a *c + b* e)

w0 <- 1
w1 <- - (a*c + b*e)
w2 <-  a*c * b*e - (d^2 * a * b)

(-w1 + c(-1,1) * sqrt(w1^2 - 4*w0*w2))/(2*w0)


lambdas <- (-w1 + c(1,-1) * sqrt(w1^2 - 4*w0*w2))/(2*w0)
ev_u
lambdas

# conditions
a*c - lambdas # not equal to zero
a*b*d^2 + (lambdas - a*c)^2 # not equal to zero

lambda <- lambdas[1]
v1 <- sqrt( (b*d^2)/(a*b*d^2 + (lambda-a*c)^2)  )
v2 <- sqrt((1 - a*v1^2)/b) # (lambda-a*c)/sqrt(a*b^2*d^2 + b*(lambda-a*c)^2)
v <- c(rep(v1, a), rep(v2, b))
ev_u; lambda
range(St %*% v - lambda*v)
V1 <- v

lambda <- lambdas[2]
v1 <- sqrt( (b*d^2)/(a*b*d^2 + (lambda-a*c)^2)  )
v2 <- - sqrt((1 - a*v1^2)/b) # (lambda-a*c)/sqrt(a*b^2*d^2 + b*(lambda-a*c)^2)
v <- c(rep(v1, a), rep(v2, b))
ev_u; lambda
range(St %*% v - lambda*v)
V2 <- v

V3 <- c(rep(c(1,-1)/sqrt(a), each = a/2), rep(0,b))
sum(V3^2)
lambdas <- c(lambdas, (a/2)*(mean(St[1,1:(a/2)]) - mean(St[1,a/2 + 1:(a/2)])))

B <- matrix(0,p,2)
B[1:a,1] <- B[a + 1:b,2] <- 1
G1 <- B %*% solve(t(B) %*% solve(St) %*% B) %*% t(B) %*% solve(St)
image2(G1)
range(B %*% solve(t(B) %*% solve(St) %*% B) %*% t(B) %*% solve(St) - B %*% MASS::ginv(B))

B <- cbind(B, t(sapply(1:p, function(k) ((1:d0) == ij[k,1]) +  ((1:d0) == ij[k,2])) ))
G2 <- B %*% MASS::ginv(t(B) %*% solve(St) %*% B) %*% t(B) %*% solve(St)
range(G2 - t(G2))
image2(G2- t(G2))
image2(G2)
range(B %*% MASS::ginv(t(B) %*% solve(St) %*% B) %*% t(B) %*% solve(St) - B %*% MASS::ginv(B))

St[1:a,1:a] %*% rep(1/a, a)

(St %*% (diag(p) - G2)) / (diag(p) - G2)
(St %*% (G2 - G1)) / (G2 - G1)

G1 <- matrix(0,p,p); G1[1:a, 1:a] <- 1/a
(St %*% G1[,1]) / G1[,1]

c(St[1:a,1:a] %*% (G2[1:a,1] - G1[1:a,1])) / (G2[1:a,1] - G1[1:a,1])


B <- matrix(0,p,3)
B[1:(a/2),1] <- B[a/2 + 1:(a/2),2] <- B[a + 1:b,3] <- 1

