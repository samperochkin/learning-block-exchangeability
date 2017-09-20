#.libPaths('/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')

library(parallel, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(MASS, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(matrixcalc, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(pcaPP, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
#library(data.table)

# load scenarios
#setwd("C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/old code/Simulation_study")
source(paste(getwd(),"/Tau_generatorsHD.R", sep = ""))


# scenario grid
n.scenario <- c(25,50,75,100)
scenario.grid <- expand.grid(list(n = n.scenario, Tau.num = 1:length(generator.list)))
save(scenario.grid , file = paste(getwd(), "/scenariogrid", sep = ""))
num.rep <- 500

print(scenario.grid)
print(nrow(scenario.grid))


# load functions
#sapply(list.files(pattern="[.]R$", path = "C:/Users/Samuel/Dropbox/Perreault-Duchesne-Neslehova/Additions/new code", full.names=TRUE), source)
sapply(list.files(pattern="[.]R$", path = paste(getwd(), "/functions", sep = ""), full.names=TRUE), source)

# create cluster
cl <- makeCluster(detectCores()-1)

# export libraries
clusterEvalQ(cl,library("MASS", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("matrixcalc", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("pcaPP", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("data.table", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))

# export objects
clusterExport(cl, ls())

for(i in 1:nrow(scenario.grid)){


  clusterExport(cl, "i")

  n <- scenario.grid[i,1]
  

  Tau <- generator.list[[scenario.grid[i,2]]]

  for(j in 1:10){
    j.vec <- 1:50 + (j - 1) * 50

    res  <- parSapply(cl, j.vec, function(k){

      set.seed(666*i + k)

      n <- scenario.grid[i,1]
      Tau <- generator.list[[scenario.grid[i,2]]]
      
      d <- nrow(Tau)

      X <- mvrnorm(n = n, rep(0,d), Sigma = sin(pi / 2 * Tau), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      Tau.hat <- cor.fk(X)
      Theta.hat <- computeTh(X)
      
      #### run Algorithm 1 with three different shrinkage parameters ####
      
      pathBuilder0(Tau.hat, Theta.hat, n)
    })
    save(res , file = paste("HDres", i, "_", j, sep = ""))
  }
    
  get.there <- TRUE
  save(get.there , file = paste("HDsuccess", i, sep = ""))
}
  

stopCluster(cl)

##########
# res <- pathBuilder0(Tau.hat, Theta.hat, n)
#
# i <- 17
# DD <- tcrossprod(res$D[[i]])
# TT <- ((DD %*% (Tau.hat - diag(d)) %*% DD) / (DD %*% (1 - diag(d)) %*% DD))
# diag(TT) <- 1
#
# par(mfrow = c(1,2), mar = c(1,1,1,1))
# image(t(TT[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))
# image(t(Tau[d:1,]), axes=FALSE, zlim=c(-1,1), col=colfunc(50))
##########
