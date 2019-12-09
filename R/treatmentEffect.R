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
  CovData <- cbind(CovData, matrix(1, nrow = nrow(CovData), 1))

  ParApt0 <- Results$estimate[paste0('Outcome.0.', c(OutcomeCov, 'Const'))]
  ParApt1 <- Results$estimate[paste0('Outcome.1.', c(OutcomeCov, 'Const'))]

  DistParEst <- calcPar(Results)

  # Treatment effect on the treated
  TreatedObs <- which(RegData[, SelectDep, with = F] == 1)
  SelectPar <- matrix(Results$estimate[paste0('Select.', c(SelectCov, 'Const'))], ncol = 1)

  SelectCovDataTreated <- as.matrix(RegData[TreatedObs, SelectCov, with = F])
  SelectCovDataTreated <- cbind(SelectCovDataTreated, matrix(1, nrow = nrow(SelectCovDataTreated), 1))
  MillsRatioTreated <- stats::dnorm(SelectCovDataTreated %*% SelectPar)/stats::pnorm(SelectCovDataTreated %*% SelectPar)

  SelectCovDataUnTreated <- as.matrix(RegData[-TreatedObs, SelectCov, with = F])
  SelectCovDataUnTreated <- cbind(SelectCovDataUnTreated, matrix(1, nrow = nrow(SelectCovDataUnTreated), 1))
  MillsRatioUnTreated <- stats::dnorm(SelectCovDataUnTreated %*% SelectPar)/(1-stats::pnorm(SelectCovDataUnTreated %*% SelectPar))

  # Treatment effect
  EY1.A1 <- CovData[TreatedObs, ] %*% matrix(ParApt1) + DistParEst['Outcome.1.Sigma', 1]*DistParEst['Outcome.1.Rho', 1]*MillsRatioTreated
  EY0.A0 <- CovData[-TreatedObs, ] %*% matrix(ParApt0) - DistParEst['Outcome.0.Sigma', 1]*DistParEst['Outcome.0.Rho', 1]*MillsRatioUnTreated

  EY0.A1 <- CovData[TreatedObs, ] %*% matrix(ParApt0) + DistParEst['Outcome.0.Sigma', 1]*DistParEst['Outcome.0.Rho', 1]*MillsRatioTreated
  EY1.A0 <- CovData[-TreatedObs, ] %*% matrix(ParApt1) - DistParEst['Outcome.1.Sigma', 1]*DistParEst['Outcome.1.Rho', 1]*MillsRatioUnTreated

  effectM <- data.frame(matrix(c(mean(EY1.A1), mean(EY0.A1),
                                 mean(EY1.A0), mean(EY0.A0),
                                 NA, NA),
                        nrow = 3, byrow = T))
  colnames(effectM) <- c('Y1 (Choose to adopt)', 'Y0 (Choose not to adopt)')
  row.names(effectM) <- c('A1 (Adopted)', 'A0 (Not adopted)', 'HeterogeneityEffect')
  effectM[3, ] <- effectM[1, ] - effectM[2, ]
  effectM$TreatEffect <- effectM[, 1] - effectM[, 2]
  effectM$Label <- c('ATT', 'ATU', 'ATH')
  effectM
}




