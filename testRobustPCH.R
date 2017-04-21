setwd("K:/FORSK-Projekt/Projekter/Scientific Projects/141_Jorne_phd/Scripts/RobustStatistics")

source("./penalizedPCH.R")
source("./simulatePCH.R")

contRate <- 0.05
n <- 2000
contN <- floor(n*contRate)
cleanN <- n - contN
  
library(pch)
library(MASS)
evalTimes <- seq(0, 250, 5)

nCov <- 4
  
X <- mvrnorm(n = cleanN, mu = rep(0, nCov), diag(nCov))

time <- rexp(cleanN, rate = 1 / exp(5 +  3 * X[,1]))
censoring <- rexp(cleanN, rate = 1 / exp(6))
timeVar <- pmin(time, censoring)
d01 <- as.numeric(time < censoring)
simData <- data.frame(timeVar, d01, X)
    
X <- mvrnorm(n = cleanN, mu = rep(0, nCov), 25 * diag(nCov))
contTime <- rexp(contN, rate = 1 / exp(7 - 2 * X[, 1]))
contCensoring <- rexp(contN, rate = 1 / exp(7))
contTimeVar <- pmin(contTime, contCensoring)
contd01 <- as.numeric(contTime < contCensoring)
contData <- data.frame(timeVar = contTimeVar, d01 = contd01, X)

simData <- rbind(simData, contData)
    
library(pec)
## find some optimal parameters
censoringFit <- survfit(Surv(timeVar, 1 - d01) ~ 1)
censoringProb <- predictSurvProb(censoringFit, newdata = simData, times =  evalTimes)
## pretty inneficient way of doing this but ok ...
censoringProb2 <- diag(predictSurvProb(censoringFit, newdata = simData, 
                                       times =  simData$timeVar))[rank(simData$timeVar)]
censoringProb2[censoringProb2 == 0] <- min(censoringProb2[censoringProb2 != 0])

## censoring using a robust estimator
censoringRob <- TML.censored(log(timeVar) ~ 1, delta = (1 - simData$d01), simData, errors = "logWeibull")
estHaz <- - censoringRob$th1 / censoringRob$v1
censoringProbRob <- matrix(exp(-(evalTimes) * exp(estHaz)), 
                            nrow = nrow(simData), ncol = length(evalTimes))
censoringProb2Rob <- exp(-(simData$timeVar) * exp(estHaz))

penalizedPathBrier(timeVar = simData$timeVar, 
                   d01 = simData$d01, 
                   X = as.matrix(simData[, ! names(simData) %in% c("timeVar", "d01")]), 
                   evalTimes = evalTimes, 
                   cProbEvalTimes =  censoringProb, 
                   cProbTimeVar = censoringProb2,
                   breaks = c(0, 75, Inf),
                   lambdaSeq = 10^(- seq(2, 4, 0.1)))

penalizedPathBrier(timeVar = simData$timeVar, 
                   d01 = simData$d01, 
                   X = as.matrix(simData[, ! names(simData) %in% c("timeVar", "d01")]), 
                   evalTimes = evalTimes, 
                   cProbEvalTimes =  censoringProbRob, 
                   cProbTimeVar = censoringProb2Rob,
                   breaks = c(0, 75, Inf),
                   lambdaSeq = 10^(- seq(2, 4, 0.1)))

## test with PCH data
X <- mvrnorm(n = cleanN, mu = rep(0, nCov), diag(nCov))

time <- simulatePCH(cleanN, X, c(-5, -2, 0,0,0, -3, 1, 0,0,0), breaks = c(0, 150, Inf))
censoring <- rexp(cleanN, rate = 1 / exp(6))
timeVar <- pmin(time, censoring)
d01 <- as.numeric(time < censoring)
simData <- data.frame(timeVar, d01, X)
breaks <- c(seq(0, 200, 25), Inf)

## penalized MLE
maxLam <- lambdaMaxMLE(timeVar = simData$timeVar, d01 = simData$d01, 
                       breaks = breaks,
                       X = as.matrix(simData[, !names(simData) %in% c("timeVar", "d01")]))

mlePath <- penalizedPathMLE(timeVar = simData$timeVar, d01 = simData$d01, 
                            breaks = breaks, 
                            X = as.matrix(simData[, !names(simData) %in% c("timeVar", "d01")]),
                            lambdaSeq = 10^seq(maxLam^(1/10),-10, by = -0.2))

cvSolution <- cvMLE(timeVar = simData$timeVar, d01 = simData$d01, 
                    breaks = breaks, 
                    X = as.matrix(simData[, !names(simData) %in% c("timeVar", "d01")]))


makePlot <- function(parameters, nCov, breaks){
  for(i in 1:(nCov + 1)){
    originalParam <-  as.numeric(paramTransMat(nIntervals = length(breaks) - 1,
                                    nParameters = nCov + 1) %*% mlePath[55, ])
    
    parameterIndices <- i + 0:(length(breaks) - 2) * (nCov + 1)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 300), ylim = c(-6, 1))
    segments(x0 = 0, y0 = originalParam[parameterIndices[1]], x1 = breaks[2])
    for(j in 2:(length(breaks) - 1)){
      segments(x0 = breaks[j], y0 = originalParam[parameterIndices[j]], x1 = breaks[j + 1])
    }
  }
}

makePlot(cvSolution, nCov, breaks)
## penalized brier
## find some optimal parameters
censoringFit <- survfit(Surv(timeVar, 1 - d01) ~ 1)
censoringProb <- predictSurvProb(censoringFit, newdata = simData, times =  evalTimes)
## pretty inneficient way of doing this but ok ...
censoringProb2 <- diag(predictSurvProb(censoringFit, newdata = simData, 
                                       times =  simData$timeVar))[rank(simData$timeVar)]
censoringProb2[censoringProb2 == 0] <- min(censoringProb2[censoringProb2 != 0])

penalizedPathBrier(timeVar = simData$timeVar, d01 = simData$d01, 
                   breaks = c(0, 150, Inf), X = as.matrix(simData[, !names(simData) %in% c("timeVar", "d01")]),
                   lambdaSeq = 10^seq(0,-10, by = -0.2), evalTimes = evalTimes, 
                   cProbEvalTimes = censoringProb,
                   cProbTimeVar = censoringProb2)

## see what happens when the data is contaminated
X <- mvrnorm(n = cleanN, mu = rep(0, nCov), 25 * diag(nCov))
contTime <- rexp(contN, rate = 1 / exp(7 - 2 * X[, 1]))
contCensoring <- rexp(contN, rate = 1 / exp(7))
contTimeVar <- pmin(contTime, contCensoring)
contd01 <- as.numeric(contTime < contCensoring)
contData <- data.frame(timeVar = contTimeVar, d01 = contd01, X)

simData <- rbind(simData, contData)