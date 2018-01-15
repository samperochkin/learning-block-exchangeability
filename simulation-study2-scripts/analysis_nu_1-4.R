# Conclusion: we need to change only values for nu2 and Table D3.

load("simulation-study2-scripts/perf1_mats.R")

mat <- t(sapply(1:60, function(i){
  j <- ceiling(i/5)
  print(j)
  bestSE.mat[i,]/hatSE.mat[j,]
}))

nu2 <- (1 - matrix(rowMeans(mat), nrow = 5))

# for main paper
nu2[4,]
#[1] 0.6521659 0.6675213 0.6732653 0.6965260 0.7951041 0.7964953 0.5071977 0.5934955 0.6311683 0.5075241 0.5186562 0.5643605

# for appendix
nu2 <- (1 - matrix(rowMeans(mat), ncol = 5, byrow = TRUE))
nu2

library(xtable)
xtable(nu2)
#\begin{table}[ht]
#\centering
#\begin{tabular}{rrrrrr}
#\hline
#& 1 & 2 & 3 & 4 & 5 \\ 
#\hline
#1 & 0.00 & 0.65 & 0.65 & 0.65 & 0.65 \\ 
#2 & 0.01 & 0.67 & 0.67 & 0.67 & 0.67 \\ 
#3 & 0.67 & 0.67 & 0.67 & 0.67 & 0.67 \\ 
#4 & 0.00 & 0.65 & 0.69 & 0.70 & 0.70 \\ 
#5 & 0.27 & 0.79 & 0.80 & 0.80 & 0.79 \\ 
#6 & 0.80 & 0.80 & 0.80 & 0.80 & 0.80 \\ 
#7 & 0.00 & 0.48 & 0.50 & 0.51 & 0.49 \\ 
#8 & 0.08 & 0.58 & 0.59 & 0.59 & 0.58 \\
#9 & 0.63 & 0.63 & 0.63 & 0.63 & 0.63 \\ 
#10 & 0.00 & 0.44 & 0.49 & 0.51 & 0.52 \\ 
#11 & 0.10 & 0.48 & 0.51 & 0.52 & 0.51 \\ 
#12 & 0.49 & 0.53 & 0.56 & 0.56 & 0.54 \\ 
#\hline
#\end{tabular}
#\end{table}

# NEED TO CLARIFY VALUE IN TABLE 1


# For Table D2 nu1
nu1 <- matrix(rowMeans(getTrue.mat), ncol = 5, byrow = TRUE)
nu1

full_table <- cbind(nu1, nu2)
xtable(full_table)
#\begin{table}[ht]
#\centering
#\begin{tabular}{rrrrrrrrrrr}
#\hline
#& 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 \\ 
#\hline
#1 & 0.00 & 0.94 & 0.95 & 0.94 & 0.92 & 0.00 & 0.65 & 0.65 & 0.65 & 0.65 \\ 
#2 & 0.02 & 1.00 & 1.00 & 1.00 & 0.99 & 0.01 & 0.67 & 0.67 & 0.67 & 0.67 \\ 
#3 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 0.67 & 0.67 & 0.67 & 0.67 & 0.67 \\ 
#4 & 0.00 & 0.41 & 0.50 & 0.52 & 0.50 & 0.00 & 0.65 & 0.69 & 0.70 & 0.70 \\ 
#5 & 0.20 & 0.94 & 0.95 & 0.95 & 0.95 & 0.27 & 0.79 & 0.80 & 0.80 & 0.79 \\ 
#6 & 1.00 & 1.00 & 1.00 & 1.00 & 1.00 & 0.80 & 0.80 & 0.80 & 0.80 & 0.80 \\ 
#7 & 0.00 & 0.21 & 0.26 & 0.29 & 0.23 & 0.00 & 0.48 & 0.50 & 0.51 & 0.49 \\ 
#8 & 0.07 & 0.76 & 0.81 & 0.81 & 0.75 & 0.08 & 0.58 & 0.59 & 0.59 & 0.58 \\ 
#9 & 0.98 & 0.99 & 0.99 & 0.99 & 0.98 & 0.63 & 0.63 & 0.63 & 0.63 & 0.63 \\ 
#10 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.00 & 0.44 & 0.49 & 0.51 & 0.52 \\ 
#11 & 0.00 & 0.00 & 0.01 & 0.00 & 0.00 & 0.10 & 0.48 & 0.51 & 0.52 & 0.51 \\ 
#12 & 0.03 & 0.08 & 0.10 & 0.11 & 0.08 & 0.49 & 0.53 & 0.56 & 0.56 & 0.54 \\ 
#\hline
#\end{tabular}
#\end{table}




# For Table D3
mat0 <- sapply(1:12, function(i){
  if(i <= 6){
    true.model <- 17
  }else{
    true.model <- 13
  }
  alphas <- readRDS(paste0("simulation-study2-results/alphas/alpha0_",i))
  
  vec <- getTrue.mat[(i-1) * 5 + 4,]
  rowSums((sapply(1:ncol(alphas), function(k){
    sapply(c(.05,.1,.25,.5), function(aa){
        rev(which(alphas[,k] > aa))[1]
    })
  }) == true.model)*vec)/sum(vec)
})

mat75 <- sapply(1:12, function(i){
  if(i <= 6){
    true.model <- 17
  }else{
    true.model <- 13
  }
  alphas <- readRDS(paste0("simulation-study2-results/alphas/alpha75_",i))
  
  vec <- getTrue.mat[(i-1) * 5 + 4,]
  rowSums((sapply(1:ncol(alphas), function(k){
    sapply(c(.05,.1,.25,.5), function(aa){
      rev(which(alphas[,k] > aa))[1]
    })
  }) == true.model)*vec)/sum(vec)
})

mat0/mat75

# We should change mat0 too.
xtable(mat0[,1:9])
#\begin{table}[ht]
#\centering
#\begin{tabular}{rrrrrrrrrr}
#\hline
#& 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 \\ 
#\hline
#1 & 0.93 & 0.96 & 0.96 & 0.61 & 0.94 & 0.95 & 0.09 & 0.51 & 0.96 \\ 
#2 & 0.86 & 0.90 & 0.92 & 0.70 & 0.91 & 0.90 & 0.19 & 0.59 & 0.92 \\ 
#3 & 0.78 & 0.77 & 0.77 & 0.72 & 0.79 & 0.76 & 0.24 & 0.67 & 0.76 \\ 
#4 & 0.54 & 0.53 & 0.50 & 0.45 & 0.50 & 0.51 & 0.44 & 0.58 & 0.58 \\ 
#\hline
#\end{tabular}
#\end{table}




load("simulation-study2-scripts/xi_mat.R")

mat <- matrix(rowMeans(xi.mat),4,12)
xtable(mat)
#\begin{table}[ht]
#\centering
#\begin{tabular}{rrrrrrrrrrrrr}
#\hline
#& 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 \\ 
#\hline
#1 & 0.61 & 0.66 & 0.67 & 0.61 & 0.78 & 0.79 & 0.35 & 0.47 & 0.61 & 0.37 & 0.39 & 0.44 \\ 
#2 & 0.61 & 0.65 & 0.66 & 0.63 & 0.78 & 0.78 & 0.38 & 0.50 & 0.61 & 0.39 & 0.41 & 0.46 \\ 
#3 & 0.59 & 0.61 & 0.63 & 0.64 & 0.75 & 0.75 & 0.41 & 0.52 & 0.59 & 0.41 & 0.42 & 0.47 \\ 
#4 & 0.53 & 0.55 & 0.54 & 0.61 & 0.70 & 0.69 & 0.43 & 0.51 & 0.53 & 0.41 & 0.42 & 0.45 \\ 
#\hline
#\end{tabular}
#\end{table}