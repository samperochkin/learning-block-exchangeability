# stopping criterion -- on train directly

library(parallel)
library(MASS)
library(matrixcalc)
library(Matrix)
library(pcaPP)

# create filepath for results
Results.filepath <- file.path("C:","Users","Samuel","Dropbox","HD-sim")

source(paste(getwd(),"/Tau_generators.R", sep = ""))
# scenario grid
n.scenario <- c(25,50,75,100)
scenario.grid <- expand.grid(list(n = n.scenario, Tau.num = 1:length(generator.list)))
num.rep <- 500

cl <- makeCluster(detectCores()-1)

# export libraries
clusterEvalQ(cl,library("MASS"))
clusterEvalQ(cl,library("matrixcalc"))
clusterEvalQ(cl,library("pcaPP"))

clusterExport(cl, ls())

the.grid <- expand.grid(1:5,1:nrow(scenario.grid))
ThMSE <- NULL

for(k in 1:dim(the.grid)[1]){
  
  j <- the.grid[k,1]
  i <- the.grid[k,2]
  
  print(c(i,j))
  
  load(paste("C:/Users/Samuel/Dropbox/HD-sim/res", i, "_", j, sep = ""))
  
  num.rep <- dim(res)[2]
  
  clusterExport(cl, "i")
  clusterExport(cl, "j")
  clusterExport(cl, "res")
  clusterExport(cl, "num.rep")
  

  ThMSE <- c(ThMSE,parSapply(cl, 1:100, function(r){
    set.seed(666*i + (j - 1) * 100 + r)
    n <- scenario.grid[i,1]
    Tau <- generator.list[[scenario.grid[i,2]]]
    d <- nrow(Tau)
    X <- mvrnorm(n = n, rep(0,d), Sigma = sin(pi / 2 * Tau), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    Tau.hat <- cor.fk(X)
    
    sum((Tau.hat - Tau)^2) / (d * (d-1) / 2)
  }))
}  

stopCluster(cl)


save(ThMSE , file = paste("C:/Users/Samuel/Dropbox/HD-sim/ThMSE", sep = ""))
