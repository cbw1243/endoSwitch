#' Endogenous switching regression
#'
#' This function calculates average treatment effects from an estimated endogenous switching regression model,
#' including average treatment effects, average treatment effects on the treated, average treatment effects on
#' the untreated.
#'
#' @param Results Estimated endogenous switching regression model
#' @param RegData an data frame. Data for running the regression analysis.
#' @param OutcomeDep character. Dependent variable in the outcome equation.
#' @param SelectDep character. Dependent variable in the Selection model. The iable must be binary (0 or 1).
#' @param OutcomeCov character vector. Independent iables in the outcome equation.
#' @param SelectCov character vector. Independent variables in the selection equation.

#'
#'@references Abdulai (2016) Impact of conservation agriculture technology on
#'household welfare in Zambia. \emph{Agricultural Economics} 47:729-741
#' (\href{https://onlinelibrary.wiley.com/doi/abs/10.1111/agec.12269}{link here})
#'
#' @return A matrix that reports the calculated average treatment effects.


treatmentEffect <- function(Results, RegData, OutcomeDep, SelectDep, OutcomeCov, SelectCov){
  # Estimate treatment effects.

  # 1. Treatment effects
  CovData <- as.matrix(RegData[, OutcomeCov, with = F])

  ParApt0 <- Results$estimate[paste0('Outcome.0.', OutcomeCov)]
  ParApt1 <- Results$estimate[paste0('Outcome.1.', OutcomeCov)]

  DistParEst <- calcPar(Results)

  # Treatment effect
  TE <- CovData %*% matrix(ParApt1 - ParApt0)
  # Average treatment effect
  ATE <- mean(TE)
  ATE_SD <- stats::sd(TE)

  # Treatment effect on the treated
  TreatedObs <- which(RegData[, SelectDep, with = F] == 1)
  SelectPar <- matrix(Results$estimate[paste0('Select.', SelectCov)], ncol = 1)

  SelectCovDataTreated <- as.matrix(RegData[TreatedObs, SelectCov, with = F])
  MillsRatioTreated <- stats::dnorm(SelectCovDataTreated %*% SelectPar)/stats::pnorm(SelectCovDataTreated %*% SelectPar)

  SelectCovDataUnTreated <- as.matrix(RegData[-TreatedObs, SelectCov, with = F])
  MillsRatioUnTreated <- -stats::dnorm(SelectCovDataUnTreated %*% SelectPar)/(1-stats::pnorm(SelectCovDataUnTreated %*% SelectPar))

  TT <- CovData[TreatedObs,] %*% matrix(ParApt1 - ParApt0) +
    (DistParEst[4, 1]*DistParEst[2, 1] - DistParEst[3, 1]*DistParEst[1, 1])*MillsRatioTreated

  ATT <- mean(TT)
  ATT_SD <- stats::sd(TT)

  # Treatment effect on the untreated
  TU <- CovData[-TreatedObs,] %*% matrix(ParApt1 - ParApt0) +
    (DistParEst[4, 1]*DistParEst[2, 1] - DistParEst[3, 1]*DistParEst[1, 1])*MillsRatioUnTreated

  ATU <- mean(TU)
  ATU_SD <- stats::sd(TU)

  effectM <-
    matrix(c(ATE, ATE_SD, ATE/ATE_SD,
             ATT, ATT_SD, ATT/ATT_SD,
             ATU, ATU_SD, ATU/ ATU_SD), nrow = 3)

  row.names(effectM) <- c('Avg treatment effect', 'Avg treatment effect on treated', 'Avg treatment effect on untreated')
  colnames(effectM) <- c('Estimate', 'Std', 'z')

  effectM
}
