contData <- readRDS(file = paste("//ithfil08/Projekt/Projekter/",
                    "Scientific Projects/141_Jorne_phd/",
                    "GeneratedData/dlbclDKContAlbumin.rds",
                    sep = ""))
source(paste("//ithfil08/Projekt/Projekter/",
             "Scientific Projects/141_Jorne_phd/",
             "Scripts/RobustStatistics/robcox.r",
             sep = ""))

contData <- contData[complete.cases(contData[, c("Albumin", "Age", "LDH")]), ]

## evaluation times
evalTimes = seq(0, 3000, 50)
## calculate the censoring probabilities
censoringFit <- survfit(Surv(timeVar, 1 - d01) ~ 1, data = contData)
censoringProb <- predictSurvProb(censoringFit, newdata = contData, times =  evalTimes)
## pretty inneficient way of doing this but ok ...
censoringProb2 <- diag(predictSurvProb(censoringFit, 
                                       newdata = contData, 
                                       times =  contData$timeVar))[rank(contData$timeVar)]
censoringProb2[censoringProb2 == 0] <- min(censoringProb2[censoringProb2 != 0])


## PCH models
breaks = c(seq(0, 1500, 150), Inf)
penalizedBrier(timeVar = contData$timeVar,
               d01 = contData$d01,
               X = as.matrix(contData[c("Albumin", "Age", "LDH")]),
               evalTimes = seq(0, 3000, 50),
               breaks = c(seq(0, 1500, 150), Inf),
               cProbEvalTimes = censoringProb,
               cProbTimeVar = censoringProb2,
               c(log(sum(contData$d01) / sum(contData$timeVar)), 0, 0, 0, rep(0, 4*(length(breaks)-2))),
               lambda = 0)

penalizedMLE(timeVar = contData$timeVar,
               d01 = contData$d01,
               X = as.matrix(contData[c("Albumin", "Age", "LDH")]),
               breaks = c(seq(0, 1500, 150), Inf),
               c(log(sum(contData$d01) / sum(contData$timeVar)), 0, 0, 0, rep(0, 4*(length(breaks)-2))),
               lambda = 0)

trimmedSolution <- trimmedCox(contData$timeVar, 
                              contData$d01, 
                              x = as.matrix(contData[, c("Albumin", "Age", "LDH")]))

which(contData$Albumin < 10)
sum(trimmedSolution$outliers %in% which(contData$Albumin < 10))


library(RobustAFT)
TML.censored(log(timeVar) ~ Age + Albumin + LDH, data = contData, contData$d01, errors = "logWeibull" )
