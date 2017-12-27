#####################

# You don't need to run this file. Everything it produces is in the repo.

#####################

library(rvest)
library(quantmod)
library(fGarch)
library(pcaPP)
library(data.table)

# Say we are interested in NASDAQ100 (in that case, we don't have to download the following ticker lists...)

#NASDAQ <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nasdaq&render=download")[,1]
#NYSE <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nyse&render=download")[,1]
#AMEX <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=amex&render=download")[,1]



# rvest the tickers of NASDAQ 100 stocks ----------------------------------

# We rvest the tickers from internet (www.cnbc.com/nasdaq-100) -- might need to be changed if the site is modified
#symbols <- read_html("https://www.cnbc.com/nasdaq-100/")

#symbols <- symbols %>%
#    html_nodes(".text a") %>%
#    html_text()

#saveRDS(symbols, "application-NASDAQ100/symbols")

symbols <- readRDS("application-NASDAQ100/symbols")



# Put necessary prefix to the tickers -------------------------------------

# In this case, the next part is a bit useless as all stocks are NASDAQ...
# But it can be used to put the right prefix to a given stock.
# Failing to do so sometimes triggers an error in
# quantmod::getSymbols.google. That happens when you hit
# a stock with a non-unique ticker (across the broader market)

#symbols.vec <- as.vector(sapply(symbols, function(s){
#  if(sum(s == NASDAQ) == 1){
#    return(as.character(paste0("NASDAQ:",s)))
#  }else if(sum(s == NYSE) == 1){
#    return(paste0("NYSE:",s))
#  }else if(sum(s == AMEX) == 1){
#    return(paste0("AMEX:",s))
#  }else{
#    print(paste0(s, " is not present anywhere..."))
#  }
#}))
#rm("NASDAQ", "NYSE", "AMEX")


# We simply do
symbols.vec <- paste0("NASDAQ:",symbols)




# Get stocks value at close from google finance ---------------------------

getSymbols.google(symbols.vec,
                  env = .GlobalEnv,
                  return.class = 'xts',
                  from = "2017-01-01",
                  to = "2017-09-30")



# Create data matrix ------------------------------------------------------

X <- do.call(cbind, sapply(symbols.vec, function(s){
  print(s)
  eval(parse(text = paste("`",s,"`", sep = "")))[,4]
}, simplify = FALSE))
rm(symbols.vec)

colnames(X) <- symbols


# Compute log-returns and get residual from GARCH(1,1) models -------------

X <- diff(as.matrix(X))/as.matrix(X[-dim(X)[1],])
X <- log(1+X)
par(mfrow = c(4,5))
sapply(1:ncol(X), function(i){
  plot(X[,i], type = "l")
})

X <- sapply(1:(dim(X)[2]),function(i){
  print(symbols[i])
  garchFit(formula = ~garch(1, 1), data = X[,i],
           cond.dist="QMLE", trace = FALSE, delta = 2, skew = 1, shape = 10)@residuals})


par(mfrow = c(4,5))
sapply(1:ncol(X), function(i){
  plot(X[,i], type = "l", main = symbols[i])
})

fwrite(data.table(X), "application-NASDAQ100/X_NASDAQ100")


# Little peek at the data -------------------------------------------------
n <- nrow(X)
d <- ncol(X)


Tau.hat <- cor.fk(X)
par(mar=c(0,0,0,0), mfrow = c(1,1))
colfunc <- colorRampPalette(c("darkred","lightyellow","forestgreen"))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim = c(-1,1), col=colfunc(100))




# Get full names, sectors and industries given on nasdaq.com --------------

# This is for further analyses of the results

NASDAQ <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nasdaq&render=download")

NASDAQ100 <- t(sapply(symbols, function(s){
  NASDAQ[which(NASDAQ[,1] == s),]
},simplify = FALSE))

NASDAQ100 <- rbindlist(NASDAQ100)

fwrite(NASDAQ100, "application-NASDAQ100/NASDAQ100")
