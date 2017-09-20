
library(parallel, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(MASS, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(matrixcalc, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')
library(pcaPP, lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3')


source(paste(getwd(),"/Tau_generatorsHD.R", sep = ""))
# scenario grid
n.scenario <- c(25,50,75,100)
scenario.grid <- expand.grid(list(n = n.scenario, Tau.num = 1:length(generator.list)))
num.rep <- 500

cl <- makeCluster(detectCores()-1)

# export libraries
clusterEvalQ(cl,library("MASS", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("matrixcalc", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))
clusterEvalQ(cl,library("pcaPP", lib.loc = '/users/sperreau/R/x86_64-redhat-linux-gnu-library/3.3'))

clusterExport(cl, ls())

#the.grid <- expand.grid(1:10, 1:nrow(scenario.grid))
# we just do the first half here!
the.grid <- expand.grid(1:10, (nrow(scenario.grid)/2 + 1):nrow(scenario.grid))

allMSEHD <- matrix(,100,0)

for(k in 1:dim(the.grid)[1]){
  
  j <- the.grid[k,1]
  i <- the.grid[k,2]
  
  print(c(i,j))
  
  load(paste("HDres", i, "_", j, sep = ""))
  
  num.rep <- dim(res)[2]
  
  d <- 50 + 50 * (i > 8)
  
  clusterExport(cl, "i")
  clusterExport(cl, "j")
  clusterExport(cl, "d")
  clusterExport(cl, "res")
  clusterExport(cl, "num.rep")
  

  allMSEHD <- cbind(allMSEHD, parSapply(cl, 1:50, function(r){
    set.seed(666*i + (j - 1) * 50 + r)
    Tau <- generator.list[[scenario.grid[i,2]]]
    Tau.hat <- res[,r]$Th
    
    sapply(1:d, function(s){
      DD <- tcrossprod(res[,r]$D[[s]])
      TT <- ((DD %*% (Tau.hat - diag(d)) %*% DD) / (DD %*% (1 - diag(d)) %*% DD))
      diag(TT) <- 1
      sum((TT - Tau)^2) / (d * (d-1) / 2)
    })
    }))
}  

stopCluster(cl)

save(allMSEHD , file = paste("allMSEHD2", sep = ""))
