simulatePCH <- function(n, X, parameters, breaks){
  nIntervals <- length(breaks) - 1
  nCoef <- ncol(X)
  
  coef <- matrix(parameters[- seq(from = 1, to = length(parameters), by = nCoef + 1)],
                 ncol = nCoef, byrow = T)
  intercept <- parameters[seq(from = 1, to = length(parameters), by = nCoef + 1)]
  
  rates <- exp(sapply(1:nIntervals, function(x) intercept[x] + X %*% coef[x, ]))
  simData <- sapply(1:nIntervals, function(x) rexp(n, rates[, x]))
  simData <- simData + matrix(breaks[- length(breaks)], byrow = T, nrow = n, ncol = nIntervals)
  simData[simData > matrix(breaks[- 1], byrow = T, ncol = nIntervals, nrow = n)] <- Inf
  do.call(pmin, as.data.frame(simData))
}