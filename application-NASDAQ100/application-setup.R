library(rvest) # Warning: you need to use the selectorgadget tool and make appropriate corrections if the site has undergone changes...
library(quantmod)
library(fGarch)
library(pcaPP)
library(data.table)

NASDAQ <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nasdaq&render=download")[,1]
NYSE <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nyse&render=download")[,1]
AMEX <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=amex&render=download")[,1]

# say we are interested in NASDAQ100
symbols <- read_html("https://www.cnbc.com/nasdaq-100/")
symbols <- symbols %>%
    html_nodes(".text a") %>%
    html_text()

# In this case, the next part is a bit useless as all stocks are NASDAQ...
# But it can be used to put the right prefix to a given stock.
# Failing to do so sometimes triggers an error in
# quantmod::getSymbols.google. That happens when you hit
# a stock with a non-unique ticker (across the broader market)
symbols.vec <- as.vector(sapply(symbols, function(s){
  if(sum(s == NASDAQ) == 1){
    return(as.character(paste0("NASDAQ:",s)))
  }else if(sum(s == NYSE) == 1){
    return(paste0("NYSE:",s))
  }else if(sum(s == AMEX) == 1){
    return(paste0("AMEX:",s))
  }else{
    print(paste0(s, " is not present anywhere..."))
  }
}))

rm("NASDAQ", "NYSE", "AMEX")

# get values at close for specified dates
getSymbols.google(symbols.vec,
                  env = .GlobalEnv,
                  return.class = 'xts',
                  from = "2017-01-01",
                  to = "2017-09-17")


# create matrix X
X <- do.call(cbind, sapply(symbols.vec, function(s){
  print(s)
  eval(parse(text = paste("`",s,"`", sep = "")))[,4]
}, simplify = FALSE))

rm(symbols.vec)
colnames(X) <- symbols


# compute log-returns and get residuals of a GARCH(1,1)
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


# set scalar values and take a look at Tau.hat
n <- nrow(X)
d <- ncol(X)

Tau.hat <- cor.fk(X)
par(mar=c(0,0,0,0), mfrow = c(1,1))
colfunc <- colorRampPalette(c("darkred","lightyellow","forestgreen"))
image(t(Tau.hat[d:1,]), axes=FALSE, zlim = c(-1,1), col=colfunc(100))

fwrite(data.table(X), "X_NASDAQ100")



# get full names, sectors, industry
NASDAQ <- read.csv("http://www.nasdaq.com/screening/companies-by-name.aspx?letter=0&exchange=nasdaq&render=download")

NASDAQ100 <- t(sapply(symbols, function(s){
  NASDAQ[which(NASDAQ[,1] == s),]
},simplify = FALSE))

NASDAQ100 <- rbindlist(NASDAQ100)

fwrite(NASDAQ100, "NASDAQ100")
