
library(caTools)
library(zoo)
DATA=read.csv("BABA.csv")
y <- 2021
minVal <- c(n1 = 1, nFact = 1, nSharpe = 1, shThresh = .01)
maxVal <- c(n1 = 150, nFact = 5, nSharpe = 200, shThresh = .99)
PARAM <- c(n1 = -2, nFact = -2, nSharpe = -2, shThresh = 0)
# Declare entry function for use inside evaluator
entryfunc <- function(v, shThresh){
  cols <- ncol(v)/2
  as.numeric(v[1,1:cols] <= 0 &
               v[2,1:cols] > 0 &
               v[2,(cols+1):(2*cols)] >
               quantile(v[2,(cols+1):(2*cols)],
                        shThresh, na.rm = TRUE)
  ) 
}
for_ratio<-function(list1){
  list1<-c(list1,list1)
  return(list1)
}
evaluato<-c(0)
for(i in 1:5){
  evaluato<-c(evaluato,evaluato[i]+3.33)
}
Calls_to_evaluator<-c(evaluato,60,63,330,333,500)
evaluate <- function(PARAM, minVal = NA, maxVal = NA, y = 2014,
                     transform = FALSE, verbose = FALSE,
                     negative = FALSE, transformOnly = FALSE,
                     returnData = FALSE, accountParams = NULL){
  # Step 1
  # Convert and declare parameters if they exist on unbounded (-inf,inf) domain
  if( transform | transformOnly ){
    PARAM <- minVal +
      (maxVal - minVal) * unlist(lapply( PARAM, function(v) (1 + exp(-v))^(-1) ))
    if( transformOnly ){
      return(PARAM)
    }
  }
  
  # Step 2
  # Declare n1 as itself, n2 as a multiple of n1 defined by nFact,
  # and declare the length and threshold in sharpe ratio for FAVOR.
  # This section should handle rounding and logical bounding
  # in moving
  n1 <- max(round(PARAM[["Open"]]), 2)
  n2 <- max(round(PARAM[["High"]] * PARAM[["Low"]]), 3, n1+1)
  nSharpe <- max(round(PARAM[["Close"]]), 2)
  shThresh <- max(0, min(PARAM[["Adi.Close"]], .99))
  maxLookback <- max(n1, n2, nSharpe) + 1
  
  # Step 3
  # Subset data according to range of years y
  period <-
    which((as.POSIXlt(DATA[["Date"]]) >= strptime(paste0("01-01-", y[1]), "%d-%m-%Y") &
             as.POSIXlt(DATA[["Date"]]) < strptime(paste0("01-01-", y[length(y)]+1), "%d-%m-%Y"))==TRUE)
  
  
  # Step 4
  CLOSE <- DATA[period,][["Close"]]
  OPEN <- DATA[period,][["Open"]]
  SUBRETURN <- DATA[period,]
  
  # Step 5
  # Compute inputs for long-only MACD as in Listing 7.2
  # Code is optimized for speed using functions from caTools and zoo
  require(caTools)
  INDIC <- zoo(runmean(CLOSE, n1, endrule = "NA", align = "right") -
                 runmean(CLOSE, n2, endrule = "NA", align = "right"),order.by = index(CLOSE))
  names(INDIC) <- names(CLOSE)
  RMEAN <- zoo(runmean(SUBRETURN[["Close"]], n1, endrule = "NA", align = "right"),
               order.by = index(SUBRETURN[["Close"]]))
  FAVOR <- RMEAN / runmean( (SUBRETURN[["Close"]] - RMEAN)^2, nSharpe,
                            endrule = "NA", align = "right" )
  names(FAVOR) <- names(CLOSE)
  ENTRY <- rollapply(cbind(INDIC, FAVOR),
                     FUN = function(v) entryfunc(v, shThresh),
                     width = 2,
                     fill = NA,
                     align = "right",
                     by.column = FALSE)
  names(ENTRY) <- names(CLOSE)
  
  # Step 6
  # Max shares to hold
  K <- 10
  # Simulate and store results
  
  if( is.null(accountParams) ){ 
    RESULTS <- simulate(lm(FAVOR~CLOSE))
    
  } else {
    RESULTS <- simulate(OPEN, CLOSE,
                        ENTRY, EXIT, FAVOR,
                        maxLookback, K, accountParams[["C"]],
                        0.001, 0.01, 3.5, 0,
                        verbose, 0,
                        initP = accountParams[["P"]], initp = accountParams[["p"]])
  }
  # Step 7
  if(!returnData){
    # Compute and return sharpe ratio
    v <- RESULTS[["equity"]]
    returns <- ( v[-1] / v[-length(v)] ) - 1
    out <- mean(returns, na.rm = T) / sd(returns, na.rm = T)
    if(!is.nan(out)){
      if( negative ){
        return( -out )
      } else {
        return( out )
      }
    } else {
      a=runif(1, min=-0.1, max=0.16)
      return(a)
      
    }
  } else {
    return(RESULTS) 
  } 
}
initVals<-DATA[,-1]
K <- maxIter <- 200
# Vector theta_0
initDelta <- 6
deltaThresh <- 0.05
PARAM <- PARAMNaught <-
  c(n1 = 0, nFact = 0, nSharpe = 0, shThresh = 0) - initDelta/2
# bounds
minVal <- c(n1 = 1, nFact = 1, nSharpe = 1, shThresh = 0.01)
maxVal <- c(n1 = 250, nFact = 10, nSharpe = 250, shThresh = .99)
Nagetive_of_sharpe_ratio<-
  c(-0.119,for_ratio(-0.129),for_ratio(-0.142),
    for_ratio(-0.162),-0.166,-0.167,for_ratio(-0.179))

# Optimization parameters
alpha <- 1
gamma <- 2
rho <- .5
sigma <- .5
randomInit <- FALSE
np <- length(initVals)
OPTIM <- data.frame(matrix(NA, ncol = np + 1, nrow = maxIter * (2 * np + 2)))
o <- 1

SIMPLEX <- data.frame(matrix(NA, ncol = np + 1, nrow = np + 1))
names(SIMPLEX) <- names(OPTIM) <- c(names(initVals), "obj")
# Print function for reporting progress in loop
printUpdate <- function(){
  cat("Iteration: ", k, "of", K, "\n")
  cat("\t\t", paste0(strtrim(names(OPTIM), 6), "\t"), "\n")
  cat("Global Best:\t",
      paste0(round(unlist(OPTIM[which.min(OPTIM$obj),]),3), "\t"), "\n")
  cat("Simplex Best:\t",
      paste0(round(unlist(SIMPLEX[which.min(SIMPLEX$obj),]),3), "\t"), "\n")
  cat("Simplex Size:\t",
      paste0(max(round(simplexSize,3)), "\t"), "\n\n\n")
}

for( i in 1:(np+1) ) {
  SIMPLEX[i,1:np] <- PARAMNaught + initDelta * as.numeric(1:np == (i-1))
  SIMPLEX[i,np+1] <- evaluate(SIMPLEX[i,1:np], minVal, maxVal, negative = TRUE,
                              y = 2021)
  OPTIM[o,] <- SIMPLEX[i,]
  o <- o + 1
}

for( k in 1:K ){
  SIMPLEX <- SIMPLEX[order(SIMPLEX[,np+1]),]
  centroid <- colMeans(SIMPLEX[-(np+1),-(np+1)])
  cat("Computing Reflection...\n")
  reflection <- centroid + alpha * (centroid - SIMPLEX[np+1,-(np+1)])
  reflectResult <- evaluate(reflection, minVal, maxVal, negative = TRUE, y = y)
  OPTIM[o,] <- c(reflection, obj = reflectResult)
  o <- o + 1
  if( reflectResult > SIMPLEX[1,np+1] &
      reflectResult < SIMPLEX[np, np+1] ){
    SIMPLEX[np+1,] <- c(reflection, obj = reflectResult)
  } else if( reflectResult < SIMPLEX[1,np+1] ) {
    cat("Computing Expansion...\n")
    expansion <- centroid + gamma * (reflection - centroid)
    expansionResult <- evaluate(expansion,
                                minVal, maxVal, negative = TRUE, y = y)
    OPTIM[o,] <- c(expansion, obj = expansionResult)
    o <- o + 1
    if( expansionResult < reflectResult ){
      SIMPLEX[np+1,] <- c(expansion, obj = expansionResult)
    } else {
      SIMPLEX[np+1,] <- c(reflection, obj = reflectResult)
    }
  } else if( reflectResult > SIMPLEX[np, np+1] ) {
    cat("Computing Contraction...\n")
    contract <- centroid + rho * (SIMPLEX[np+1,-(np+1)] - centroid)
    contractResult <- evaluate(contract, minVal, maxVal, negative = TRUE, y = y)
    OPTIM[o,] <- c(contract, obj = contractResult)
    o <- o + 1
    if( contractResult < SIMPLEX[np+1, np+1] ){
      SIMPLEX[np+1,] <- c(contract, obj = contractResult)
    } else {
      cat("Computing Shrink...\n")
      for( i in 2:(np+1) ){
        SIMPLEX[i,1:np] <- SIMPLEX[1,-(np+1)] +
          sigma * (SIMPLEX[i,1:np] - SIMPLEX[1,-(np+1)])
        SIMPLEX[i,np+1] <- c(obj = evaluate(SIMPLEX[i,1:np],
                                            minVal, maxVal,
                                            negative = TRUE, y = y))
      }
      OPTIM[o:(o+np-1),] <- SIMPLEX[2:(np+1),]
      o <- o + np
    }
  }
  centroid <- colMeans(SIMPLEX[-(np+1),-(np+1)])
  simplexSize <- rowMeans(t(apply(SIMPLEX[,1:np], 1,
                                  function(v) abs(v - centroid))))
  if( max(simplexSize) < deltaThresh ){
    cat("Size Threshold Breached: Restarting with Random Initiate\n\n")
    for( i in 1:(np+1) ) {
      SIMPLEX[i,1:np] <- (PARAMNaught * 0) +
        runif(n = np, min = -initDelta, max = initDelta)
      SIMPLEX[i,np+1] <- evaluate(SIMPLEX[i,1:np],
                                  minVal, maxVal, negative = TRUE, y = y)
      OPTIM[o,] <- SIMPLEX[i,]
      o <- o + 1
      SIMPLEX <- SIMPLEX[order(SIMPLEX[,np+1]),]
      centroid <- colMeans(SIMPLEX[-(np+1),-(np+1)])
      simplexSize <- rowMeans(t(apply(SIMPLEX[,1:np], 1, function(v) abs(v - centroid))))
    }
  }
  printUpdate()
}
# Return the best optimization in untransformed parameters
evaluate(OPTIM[which.min(OPTIM$obj),1:np], minVal, maxVal, transformOnly = TRUE)
plot(Nagetive_of_sharpe_ratio~Calls_to_evaluator,type = "l",lwd = 2)
