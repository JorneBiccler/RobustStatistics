## load the necessary packages
library(RobustAFT)

#########################
### GENERAL FUNCTIONS ###
#########################


firstIndexIncInterval <- function(evalTimes, breaks){
  ## function that returns the last index of the breaks vector
  ## for which the evalTimes value is larger than the corresponding breaks value
  tempMatrix <- matrix(rep(breaks, length(evalTimes)), ncol = length(breaks), byrow = T)
  rowSums(tempMatrix <= evalTimes)
}

predictPCH <- function(parameters, breaks, X, evalTimes){
  ## Function to calculate the predicted survival of new observations
  ## at the times provided.
  ## ARGS:
  ##   parameters: a vector containing the intercepts and coefficients,
  ##               the structure should be as follows: first the first intercept
  ##               then the coefficients corresponding to the first interval
  ##               after which the same is done for the second interval etc.
  ##   breaks: the breakpoints 
  ##   X: covariate matrix
  ##   evalTimes: times for which predictions should be returned
  
  nCoef <- ncol(X)
  nIntervals <- length(breaks) - 1
  
  ## split the parameters into a matrix of coefficients, with columns corresponding
  ## to coefficients and rows to intervals, and a vector of intercepts
  coef <- matrix(parameters[- seq(from = 1, to = length(parameters), by = nCoef + 1)],
                 ncol = nCoef, byrow = T)
  intercept <- parameters[seq(from = 1, to = length(parameters), by = nCoef + 1)]
  
  ## case of an exponential model (i.e. breaks = c(0, Inf) or length(breaks) == 2)
  if(length(breaks) == 2){
    estCumHaz <- do.call(cbind, lapply(1:length(evalTimes),
                                       function(x) (evalTimes[x] - breaks[1]) *
                                         exp(intercept[1] +
                                               as.numeric(X %*% coef[1, ]))))
    
    ## jump out of the function and return the survival estimates
    predicted <- exp(- estCumHaz)
    return(predicted)
  }
  
  ## fill a matrix with the hazard contribution of the intervalls that are fully
  ## covered at the evaluation times
  fullHaz <- matrix(0, 
                    nrow = nrow(X), 
                    ncol = (length(breaks) - 2))
  
  for(j in 1:(length(breaks) - 2)){
    fullHaz[, j] <- as.numeric((breaks[j + 1] - breaks[j]) * 
                                 exp(intercept[j] + X %*% coef[j, ]))
  }
  
  ## sum the different chunks into a cumulative hazards
  fullHaz <- t(apply(fullHaz, 1, cumsum))
  if(length(breaks) == 3){
    ## if there are only 2 intervals, the result of the apply is a vector,
    ## to make everything consistent this is turned into the correct matrix
    fullHaz <- t(fullHaz)
  }
  
  ## create some help-vectors to show to which interval the different evaluation times belong
  intPlus1 <- firstIndexIncInterval(evalTimes, breaks)
  zeroIntervalsCovered <- sum(intPlus1 == 1) 
  fullyCoveredIntervals <- intPlus1[!(intPlus1 == 1)]
  
  ## calculate the cumulative hazard
  estCumHaz <- cbind(matrix(0, 
                            nrow = nrow(X),
                            ncol = zeroIntervalsCovered), 
                     fullHaz[, fullyCoveredIntervals - 1]) +
    do.call(cbind, lapply(1:length(evalTimes),
                          function(x) (evalTimes[x] - breaks[intPlus1[x]]) * 
                            exp(intercept[intPlus1[x]] +
                                  as.numeric(X %*% coef[intPlus1[x], ]))))
  
  ## return the survival estimates
  predicted <- exp(-estCumHaz)
  predicted
}

timeAtRisk <- function(evalTimes, breaks, ncolRep, nrowRep){
  ## function that creates a matrix with as components the time in the interval.
  ## the rows correspond with the evalTimes the columns with the intervals
  ## each row is repeated nrowRep times (usually the number of observations)
  ## and each column is repeated ncolRepTimes (usually the number of covariates)
  
  diffMatrix <- matrix(rep(diff(breaks), length(evalTimes)), byrow = T, nrow = length(evalTimes))
  colIndices <- firstIndexIncInterval(evalTimes, breaks)
  
  ## there might be a more efficient (vectorized) way to do this:
  diffMatrix[cbind(1:length(evalTimes), colIndices)] <- evalTimes - breaks[colIndices]
  for(i in 1:length(evalTimes)){
    if(colIndices[i] + 1 < length(breaks)){
      diffMatrix[i, (colIndices[i] + 1):ncol(diffMatrix)] <- 0
    }
  }
  
  diffMatrix <- diffMatrix[rep(1:nrow(diffMatrix), each = nrowRep), ]
  diffMatrix[, rep(1:ncol(diffMatrix), each = ncolRep)]
}

extendedY <- function(evalTimes, timeVar){
  ## function that transforms the time variable into 
  ## a new outcome variable with length(evaltimes) * lenght(timeVar) entries
  ## such that the entries are as specified in the least squares equations.
  
  rep(timeVar, length(evalTimes)) > rep(evalTimes, each = length(timeVar))
}

paramTransMat <- function(nIntervals, nParameters){
  ## returns a matrix to go from the differences parameterization to
  ## the standard parameterization
  matrixOfOnes <- kronecker(matrix(1, ncol = nIntervals, nrow = nIntervals),
                            diag(nParameters))
  matrixOfOnes[upper.tri(matrixOfOnes)] <- 0
  matrixOfOnes
}

getBrierScore <- function(timeVar, d01, X, cProbEvalTimes, cProbTimeVar, evalTimes, 
                          parameters, breaks, difParam = T){
  ## A function that returns the Brier score for a PCH model
  ## the censoringProb should include estimations of the censoring probability with
  ## rows corresponding to observations and columns to the evalTimes
  if(difParam){
    parameters <-  paramTransMat(nIntervals = length(breaks) - 1,
                                 nParameters = ncol(X) + 1) %*% parameters
  }
  ## calculate an "extended" Y to calculate quadratic errors from
  exY <- extendedY(evalTimes, timeVar)
  ## get the survival estimates given by a pch model at the evalTimes
  predictedSurvival <- predictPCH(parameters, breaks, X, evalTimes)
  
  ## calculate the IPW weights
  IPWWeights <- (rep(timeVar, length(evalTimes)) <= rep(evalTimes, each = length(timeVar)))^
    (1 - rep(d01, length(evalTimes)))
  IPWWeights <- IPWWeights / pmax(as.numeric(cProbEvalTimes), cProbTimeVar)
  
  ## us the quadratic representation of hte Brier score to just get a mean of squared errors
  mean(IPWWeights * (exY - predictedSurvival)^2)
}

getPenalty <- function(param, nParam){
  ## function that calculates the penalty function
  ## the parameters should be provided in the differences
  ## parameterization
  
  nInterval <- length(param) / nParam
  penalty <- 0
  for(i in 1:nParam){
    paramIndices <- i + (1:(nInterval - 1)) * nParam
    penalty <- penalty + sqrt(sum(param[paramIndices]^2))
  }
  penalty
}

proximalStep <- function(curSol, nParam, gradient, stepsize, lambda){
  ## Perform a proximal gradient descent step corresponding to a group lasso of differences penalty
  ## the curSol is the current solution (in differences representation),
  ## the gradient is the gradient (in differences representation),
  ## the stepsize is the stepzise that should be taken during the proximal gradient descent step
  ## lambda is the penalization parameter
  
  nInterval <- length(curSol) / nParam
  newSol <- curSol - stepsize * gradient
  
  if(lambda != 0){
    for(i in 1:nParam){
      paramIndices <- i + (1:(nInterval - 1)) * nParam
      stepNorm <- sqrt(sum(newSol[paramIndices]^2))
      if(stepNorm <= lambda * stepsize){
        newSol[paramIndices] <- 0
      } else {
        newSol[paramIndices] <- newSol[paramIndices] - 
          lambda * stepsize * (newSol[paramIndices]) / stepNorm
      }
    }
  }
  newSol
}


#############################################
### PCH MODEL WITH A PENALIZED BRIER LOSS ###
#############################################


brierGradient <- function(IPWWeights, timeAtRiskMat, exY, evalTimes, transformMat, 
                          largeX, X, X2, parameters, breaks){
  ## returns the gradient of the Brier score
  
  ## transform the values to the original parameterizaiton
  oPar <- transformMat %*% parameters
  ## get the predicted survival at the evaluation times
  predictedSurvival <- as.numeric(predictPCH(oPar, breaks, X, evalTimes))

  ## get the gradient of the original parameterization
  gradient <- NULL
  for(i in 1:(length(breaks) - 1)){
  gradient <- c(gradient, t(IPWWeights * (exY - predictedSurvival) * predictedSurvival * timeAtRiskMat[, i] *
                  as.numeric(exp(X2 %*% oPar[1:(ncol(X2)) + (i - 1) * ncol(X2)]))) %*% largeX)
  }
  
  ## transform the gradient to the gradient of the differences parameterization
  ## and return the vector
  t(transformMat) %*% gradient / (length(IPWWeights))
}

penalizedBrier <- function(timeVar, d01, X, evalTimes, 
                           cProbEvalTimes, cProbTimeVar,
                           breaks, init, 
                           lambda, 
                           epsilon = 10^(-5),
                           stepSize = 0.01,
                           bbSteps = T,
                           maxStep = 10^(5),
                           minStep = 10^(-5)){
  ## Returns the penalized Brier solution.
  
  # init is the initial solution (in the differences parameterization)
  # perform the first step and set-up some vectors
  oldSol <- init
  initStep <- stepSize
  
  ## calculate some data that's used in the gradient
  exY <- extendedY(evalTimes = evalTimes, 
                   timeVar = timeVar)
  IPWWeights <- exY^(1 - d01)
  IPWWeights <- IPWWeights / pmax(as.numeric(cProbEvalTimes), cProbTimeVar)
  timeAtRiskMat <- timeAtRisk(evalTimes, breaks, 1, length(timeVar))

  transformMat <- paramTransMat(nIntervals = length(breaks) - 1,
                                nParameters = ncol(X) + 1)
  X2 <- cbind(1, X)
  largeX <- X2[rep(1:length(timeVar), length(evalTimes)), ]
  
  ## get the gradient of the brier score at the initial solution
  gradientOld <- brierGradient(IPWWeights = IPWWeights, 
                               timeAtRiskMat = timeAtRiskMat, 
                               exY = exY, 
                               evalTimes = evalTimes, 
                               transformMat = transformMat, 
                               largeX = largeX,
                               X = X, 
                               X2 = X2,
                               parameters = oldSol, 
                               breaks =  breaks)
  
  ## perform a the proximal step
  newSol <- proximalStep(oldSol,
                         ncol(X) + 1,
                         gradientOld,
                         stepSize,
                         lambda)
  
  gradientCur <- gradientOld
  solDif <- epsilon + 1
  
  ## keep performing the proximal gradient steps untill a convergence
  ## criterion is met
  while((max(abs(gradientCur)) > epsilon) & sum(abs(solDif)) > epsilon){
    #  Barzilai-Borwein steps
    gradientCur <- brierGradient(IPWWeights = IPWWeights, 
                                 timeAtRiskMat = timeAtRiskMat, 
                                 exY = exY, 
                                 evalTimes = evalTimes, 
                                 transformMat = transformMat, 
                                 largeX = largeX,
                                 X = X,
                                 X2 = X2, 
                                 parameters = newSol, 
                                 breaks =  breaks)
    
    solDif <- newSol - oldSol
    if(bbSteps){
      gradientDif <- gradientCur - gradientOld
      gradientOld <- gradientCur
      stepSize <- sum(gradientDif * solDif) / sum(gradientDif^2)
      if(is.na(stepSize) | stepSize < 0){
        stepSize <- initStep
      } else if(stepSize < minStep) {
        stepSize <- minStep
      } else if(stepSize > maxStep){
        stepSize <- maxStep
      }
    }
    oldSol <- newSol
    newSol <- proximalStep(newSol,
                           ncol(X) + 1,
                           gradientCur,
                           stepSize,
                           lambda)
    solDif <- newSol - oldSol
  }
  newSol
}



## way of starting a good starting value for the optimization algorithm
trimmedStarts <- function(timeVar, d01, X, cProbEvalTimes, 
                          cProbTimeVar, evalTimes, ntry = 100, 
                          seed = NULL, alpha = 0.8){
  if(!is.null(seed)){
    set.seed(seed)
  }
  tempFrame <- data.frame(timeVar, d01, X)

  bestSol <- - survreg(formula = Surv(timeVar, d01) ~., tempFrame, dist = "exponential")$coef
  bestBrier <- getBrierScore(timeVar, d01, X, cProbEvalTimes, 
                             cProbTimeVar, bestSol, breaks = c(0, Inf), 
                             difParam = F, evalTimes = evalTimes)
  for(i in 1:ntry){
    curSample <- sample(length(timeVar), size = length(timeVar) * alpha)
    curSol <- - survreg(formula = Surv(timeVar, d01) ~., tempFrame[curSample, ], dist = "exponential")$coef
    curBrier <- getBrierScore(timeVar, d01, X, cProbEvalTimes, 
                              cProbTimeVar, curSol, breaks = c(0, Inf), 
                              difParam = F, evalTimes = evalTimes)
    if(curBrier < bestBrier){
      bestBrier <- curBrier
      bestSol <- curSol
    }
  } 
  bestSol
}

penalizedPathBrier <- function(timeVar, d01, X, evalTimes, 
                               cProbEvalTimes, cProbTimeVar,
                               breaks, lambdaSeq, ...){
  ## calculate a path of solutions, this helps with convergence issues
  ## by using warm-starts and a "smart" initial solution
  ## the result is a matrix in which each row contains esitmated coefficients
  ## for the corresponding lambda value

  # create a matrix in which the coefficients are saved
  solutions <- matrix(0, nrow = length(lambdaSeq), ncol = (ncol(X) + 1) * (length(breaks) - 1))

  # initSol <- trimmedStarts(timeVar = timeVar, d01 = d01, X = X, cProbEvalTimes = cProbEvalTimes,
  #                       cProbTimeVar = cProbTimeVar, evalTimes = evalTimes, ntry = 100)
  
  ## use a robust estimate of a weibull model to initialize the path
  TMLObj <- TML.censored(log(timeVar) ~ ., 
                         simData$d01, 
                         data = simData[, !names(simData) %in% "d01"], 
                         errors = "logWeibull")
  ## transform the robust estimates into the parameterization used in the
  ## rest of the functions
  initSol <- - TMLObj$th1 / TMLObj$v1
  
  ## get the first solution using the robust estimates as intial values
  solutions[1, ] <- penalizedBrier(timeVar = timeVar, d01 = d01, X = X, cProbEvalTimes = cProbEvalTimes,
                                   cProbTimeVar = cProbTimeVar, evalTimes = evalTimes, breaks = breaks,
                                   init = c(initSol, rep(0, (ncol(X) + 1)* (length(breaks) - 2))), 
                                   lambda = lambdaSeq[1])

  ## calculate the rest of the path by using warm starts
  ## i.e. use the solution of hte previous lambda value
  ## as starting value 
  for(i in 2:length(lambdaSeq)){
    solutions[i, ] <- penalizedBrier(timeVar = timeVar, d01 = d01, X = X, cProbEvalTimes = cProbEvalTimes,
                                     cProbTimeVar = cProbTimeVar, evalTimes = evalTimes, breaks = breaks,
                                     init = solutions[i - 1, ], lambda = lambdaSeq[i])
  }
  solutions
}

####################################
### PENALIZED PCH MODEL WITH MLE ###
####################################

mleGradient <- function(O, R, X2, breaks, transformMat, parameters){
  ## function that calculates the gradient of the log likelihood
  
  ## calculate the parameters in the standard parameterization
  originalParameters <- transformMat %*% parameters
  
  # get the gradient by using the original parameterization
  gradient <- NULL
  for(j in 1:(length(breaks) - 1)){
    gradient <- c(gradient,
                  t(O[, j]) %*% X2 - 
                    t(R[, j]) %*% (as.numeric(exp(X2 %*% originalParameters[1:ncol(X2) + (j - 1) * ncol(X2)])) 
                                   * X2))
  }
  
  # transform the gradient into the gradient of the differences parameterization
  - 1 * t(transformMat) %*% gradient / (nrow(X2))
}

penalizedMLE <- function(timeVar, d01, X, 
                         breaks, init, 
                         lambda, epsilon = 10^(-8),
                         stepSize = 0.01,
                         bbSteps = T){
  # init is the initial solution (in the differences parameterization)
  # perform the first step and set-up some vectors
  oldSol <- init
  
  # create some objects needed in the gradient
  R <- timeAtRisk(timeVar, breaks, 1, 1)
  lastIndex <- firstIndexIncInterval(timeVar, breaks)
  O <- matrix(0, ncol = ncol(R), nrow = nrow(R))
  O[cbind(1:length(timeVar), lastIndex)] <- d01
  X2 <- cbind(1, X)
  
  nIntervals <- length(breaks) - 1 
  nParameters <- length(parameters) / nIntervals
  transformMat <- paramTransMat(nIntervals,
                                nParameters)
  
  ## get the gradient at the provided intial solution
  gradientOld <- mleGradient(O, R, X2, breaks, transformMat, oldSol)
  
  ## perform the proximal gradient step
  newSol <- proximalStep(oldSol,
                         ncol(X) + 1,
                         gradientOld,
                         stepSize,
                         lambda)
  solDif <- epsilon + 1
  gradientCur <- epsilon + 1

  ## perform the proximal gradient steps untill the convergence
  ## criterion is met
  while((max(gradientCur^2) > epsilon) & sum(solDif^2) > epsilon){
    #  Barzilai-Borwein steps
    gradientCur <- mleGradient(O, R, X2, breaks, transformMat, newSol)
    solDif <- newSol - oldSol
    if(bbSteps){
      gradientDif <- gradientCur - gradientOld
      gradientOld <- gradientCur
      stepSize <- sum(gradientDif * solDif) / sum(gradientDif^2)
    }
    oldSol <- newSol
    newSol <- proximalStep(newSol,
                           ncol(X) + 1,
                           gradientCur,
                           stepSize,
                           lambda)
    solDif <- newSol - oldSol
  }
  newSol
}

penalizedPathMLE <- function(timeVar, d01, X, breaks, 
                             lambdaSeq){
  ## create a path of penalized MLE's.
  ## This tends to be more efficient than calculating them one by one
  ## since warm starts etc. can be used. The result is a matrix with rows
  ## corresponding to the solutions of different lambda values
  
  ## create a matrix in which the path will be saved
  solutions <- matrix(0, nrow = length(lambdaSeq), ncol = (ncol(X) + 1) * (length(breaks) - 1))
  
  ## use the standards survreg package to obtain estimates of a constant hazards (exponential) model
  ## this function is more stable and faster than using the proximal gradient method.
  tempFrame <- data.frame(timeVar, d01, X)
  initSol <- - coef(survreg(Surv(timeVar, d01) ~. , data = tempFrame, dist = "exponential"))
  
  ## use the MLE of the expenential model as initial values for the first solution
  solutions[1, ] <- penalizedMLE(timeVar = timeVar, d01 = d01, X = X, breaks = breaks,
                                 init = c(initSol, rep(0, (ncol(X) + 1)* (length(breaks) - 2))), 
                                 lambda = lambdaSeq[1])
  ## use warm-starts to calculate the solutions for the other lambda values
  for(i in 2:length(lambdaSeq)){
    solutions[i, ] <- penalizedMLE(timeVar = timeVar, d01 = d01, X = X, breaks = breaks,
                                   init = solutions[i - 1, ], lambda = lambdaSeq[i])
  }
  solutions
}


lambdaMaxMLE <- function(timeVar, d01, X, breaks){
  ## calculate the smalles lambda values for which the model
  ## has all the differences equal to zero, this is a good
  ## value from which to start a solution path.
  
  nParam <- (ncol(X) + 1)
  
  # create some objects needed in the gradient
  R <- timeAtRisk(timeVar, breaks, 1, 1)
  lastIndex <- firstIndexIncInterval(timeVar, breaks)
  O <- matrix(0, ncol = ncol(R), nrow = nrow(R))
  O[cbind(1:length(timeVar), lastIndex)] <- d01
  X2 <- cbind(1, X)
  
  nIntervals <- length(breaks) - 1 
  nParameters <- length(parameters) / nIntervals
  transformMat <- paramTransMat(nIntervals,
                                nParameters)
  ## fit a standard constant hazards model (exponential distribution).
  tempFrame <- data.frame(timeVar, d01, X)
  initSol <- - coef(survreg(Surv(timeVar, d01) ~. , data = tempFrame, dist = "exponential"))
  
  ## calculate the gradient of the MLE
  mleGrad <- mleGradient(O, R, X2, breaks, transformMat,
                             c(initSol, rep(0, (nParam)* (length(breaks) - 2))))
  ## get the lambda value
  lambdaMax <- 0
  for(i in 1:nParam){
    paramIndices <- i + (1:(nIntervals - 1)) * nParam
    gradientNorm <- sqrt(sum(mleGrad[paramIndices]^2))
    if(gradientNorm > lambdaMax){
      lambdaMax <- gradientNorm
    }
  }
  lambdaMax
}

cvMLE <- function(timeVar, d01, X, breaks, 
                  nLambda = 20, logLamdaRatio = 5, nFolds = 10, ...){
  ## funtion that selects the best lambda value by using a cross validated
  ## Brier score as evaluation criterion.
  
  ## getthe maximum lambda value of interest and create a sequence of lambda values
  logMaxLambda <- - log(lambdaMaxMLE(timeVar = timeVar, d01 = d01, breaks = breaks, X = X))
  lambdaSeq <- exp(- seq(logMaxLambda * (1 - logMaxLambda/nLambda) , logMaxLambda*logLamdaRatio, length.out = nLambda))
  
  # sample CV indices
  randSample <- sample(nrow(X))
  CVIndicesList <- split(randSample, ceiling(seq_along(randSample)/(length(randSample)/nFolds))) 
  
  # create a vector in which to save the brier scores
  brierScores <- rep(0, nLambda)
  # perform the CV steps
  for(i in 1:nFolds){
    trTimeVar <- timeVar[- CVIndicesList[[i]]]
    testTimeVar <- timeVar[CVIndicesList[[i]]]
    
    trd01 <- d01[- CVIndicesList[[i]]]
    testd01 <- d01[CVIndicesList[[i]]]
    
    trX <- X[- CVIndicesList[[i]], ]
    testX <- X[CVIndicesList[[i]], ]
  
    foldPath <- penalizedPathMLE(timeVar = trTimeVar, d01 = trd01, 
                                 X = trX, breaks = breaks,
                                 lambdaSeq = lambdaSeq)
    for(j in 1:nLambda){
      brierScores[j] <- brierScores[j] + getBrierScore(testTimeVar, 
                                                       d01 = testd01,
                                                       parameters = foldPath[j, ],
                                                       evalTimes = evalTimes,
                                                       X = testX, 
                                                       cProbEvalTimes = cProbEvalTimes[CVIndicesList[[i]]],
                                                       cProbTimeVar = cProbTimeVar[CVIndicesList[[i]]],
                                                       breaks = breaks,
                                                       difParam = T)
    }
  }
  ## calculate the solution of the best lambda value
  ## the intial parameter values are the ones obtained in the
  ## last CV step and for the corresponding lambda value
  penalizedMLE(timeVar = timeVar, 
               d01 = d01, 
               X = X, 
               breaks = breaks,
               init = foldPath[which.min(brierScores), ], 
               lambda = lambdaSeq[which.min(brierScores)])
}