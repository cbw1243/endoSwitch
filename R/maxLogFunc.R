maxLogFunc <- function(param){
  sigma0 <- exp(param[TotParNum - 3]); sigma1 <- exp(param[TotParNum - 2])
  rho0 <- (exp(2*param[TotParNum - 1]) - 1)/(exp(2*param[TotParNum - 1]) + 1)
  rho1 <- (exp(2*param[TotParNum]) - 1)/(exp(2*param[TotParNum]) + 1)

  SelParM <- matrix(param[1:SelParNum], SelParNum, 1)
  SelData <- cbind(as.matrix(RegData[, SelCovVar, with = F]), matrix(1, nrow(RegData), 1))
  SelCovSum <- SelData %*% SelParM

  SelLabel <- which(RegData[, SelDepVar, with = F] == 1)
  # Not treated
  ManParM0 <- matrix(param[(SelParNum + 1): (SelParNum + ManParNum)], ManParNum, 1)
  ManData0 <- cbind(as.matrix(RegData[-SelLabel, ManCovVar, with = F]),
                    matrix(1, nrow(RegData)-length(SelLabel), 1))
  Man0Sum <- as.vector(ManData0 %*% ManParM0)

  # Treated
  ManParM1 <- matrix(param[(SelParNum + ManParNum + 1): (SelParNum + 2*ManParNum)], ManParNum, 1)
  ManData1 <- cbind(as.matrix(RegData[SelLabel, ManCovVar, with = F]),
                    matrix(1, length(SelLabel), 1))
  Man1Sum <- as.vector(ManData1 %*% ManParM1)

  SelCovSum0 <- unlist(SelCovSum[-SelLabel, 1]) # Not Treated
  SelCovSum1 <- unlist(SelCovSum[SelLabel, 1]) # Treated

  ManRes0 <- unlist(RegData[-SelLabel, ManDepVar, with = F]) - Man0Sum
  ManRes1 <- unlist(RegData[SelLabel, ManDepVar, with = F]) - Man1Sum

  eta1 <- (SelCovSum1 + rho1*ManRes1/sigma1)/sqrt(1-rho1^2)
  eta0 <- (SelCovSum0 + rho0*ManRes0/sigma0)/sqrt(1-rho0^2)

  LogLike <- sum(log(pnorm(eta1)) + log(dnorm(ManRes1/sigma1)/sigma1)) + sum(log(1 - pnorm(eta0)) + log(dnorm(ManRes0/sigma0)/sigma0))
  if(isTRUE(verbose)) cat('Log likelihood value is:',LogLike, '\r')
  LogLike
}
