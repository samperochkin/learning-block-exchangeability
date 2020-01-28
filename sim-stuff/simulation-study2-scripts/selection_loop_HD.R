# stopping criterion -- on train directly

library(parallel, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(MASS, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(matrixcalc, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(Matrix, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(pcaPP, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')

# Note: package matrixcalc most likely not needed anymore...

source(paste(getwd(),"/Tau_generatorsHD.R", sep = ""))
# scenario grid
n.scenario <- c(25,50,75,100)
scenario.grid <- expand.grid(list(n = n.scenario, Tau.num = 1:length(generator.list)))
num.rep <- 500


#cl <- makeCluster(detectCores()-1)
cl <- makeCluster(11)

# export libraries
clusterEvalQ(cl,library("MASS", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("matrixcalc", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("Matrix", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("pcaPP", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))

# load stopping function

source("functions/stoppingCriterion0.R")

clusterExport(cl, ls())

# Again just dor the d = 50 cases for the moment
the.grid <- expand.grid(1:10, (nrow(scenario.grid)/2 + 1):nrow(scenario.grid))
selAll <- matrix(,100,0)
#selAll <- matrix(,100,0)

for(k in 1:dim(the.grid)[1]){
  
  j <- the.grid[k,1]
  i <- the.grid[k,2]
  
  print(c(i,j))
  
  load(paste("HDres", i, "_", j, sep = ""))
  
  num.rep <- dim(res)[2]
  
  clusterExport(cl, "i")
  clusterExport(cl, "j")
  clusterExport(cl, "res")
  clusterExport(cl, "num.rep")
  
  selAll <- cbind(selAll,parSapply(cl, 1:num.rep, function(r){
    Tau.hat <- res[,r]$Th
    
    stoppingCriterion0(Tau.hat, res[,r]$D, res[,r]$S)
  }))
    
  rm(res)
}  
save(selAll , file = paste("selAllHD2", sep = ""))

stopCluster(cl)
