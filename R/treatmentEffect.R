#' Endogenous switching regression
#'
#' This function calculates treatment effects from an estimated endogenous switching regression model.
#'
#' @param Results Estimated endogenous switching regression model
#' @param RegData Data for doing the regression analysis
#' @param ManDepVar Dependent variable in the main model
#' @param SelDepVar Dependent variable in the selection model. This should be binary (0 or 1).
#' @param ManCovVar Independent variables in the main model.
#' @param SelCovVar Independent variable in the selection model
#'
#' @return A list containing regression results.
#'
#' @export
#' @examples
#' data(ImpactData)
#' ManDepVar <- 'Output'
#' SelDepVar <- 'CA'
#' ManCovVar <- c('Age')
#' SelCovVar <- c('Age', 'Perception')
#' Results <- endoSwitch2Stage(ImpactData, ManDepVar, SelDepVar, ManCovVar, SelCovVar)
#' treatmentEffect(Results, ImpactData, ManDepVar, SelDepVar, ManCovVar, SelCovVar)

treatmentEffect <- function(Results, RegData, ManDepVar, SelDepVar, ManCovVar, SelCovVar){
  # Estimate treatment effects.

  # 1. Treatment effects
  CovVarData <- as.matrix(RegData[, ManCovVar, with = F])

  ParApt0 <- Results$estimate[paste0('Main.0.', ManCovVar)]
  ParApt1 <- Results$estimate[paste0('Main.1.', ManCovVar)]

  DistParEst <- calcPar(Results)

  # Treatment effect
  TE <- CovVarData %*% matrix(ParApt1 - ParApt0)
  # Average treatment effect
  ATE <- mean(TE)
  ATE_SD <- sd(TE)

  # Treatment effect on the treated
  TreatedObs <- which(ImpactData[, SelDepVar, with = F] == 1)
  SelPar <- matrix(Results$estimate[paste0('Select.', SelCovVar)], ncol = 1)

  SelCovDataTreated <- as.matrix(ImpactData[TreatedObs, SelCovVar, with = F])
  MillsRatioTreated <- dnorm(SelCovDataTreated %*% SelPar)/pnorm(SelCovDataTreated %*% SelPar)

  SelCovDataUnTreated <- as.matrix(ImpactData[-TreatedObs, SelCovVar, with = F])
  MillsRatioUnTreated <- -dnorm(SelCovDataUnTreated %*% SelPar)/(1-pnorm(SelCovDataUnTreated %*% SelPar))

  TT <- CovVarData[TreatedObs,] %*% matrix(ParApt1 - ParApt0) +
    (DistParEst[4, 1]*DistParEst[2, 1] - DistParEst[3, 1]*DistParEst[1, 1])*MillsRatioTreated

  ATT <- mean(TT)
  ATT_SD <- sd(TT)

  # Treatment effect on the untreated
  TU <- CovVarData[-TreatedObs,] %*% matrix(ParApt1 - ParApt0) +
    (DistParEst[4, 1]*DistParEst[2, 1] - DistParEst[3, 1]*DistParEst[1, 1])*MillsRatioUnTreated

  ATU <- mean(TU)
  ATU_SD <- sd(TU)

  effectM <-
    matrix(c(ATE, ATE_SD, ATE/ATE_SD,
             ATT, ATT_SD, ATT/ATT_SD,
             ATU, ATU_SD, ATU/ ATU_SD), nrow = 3)

  row.names(effectM) <- c('Avg treatment effect', 'Avg treatment effect on treated', 'Avg treatment effect on untreated')
  colnames(effectM) <- c('Estimate', 'Std', 'z')

  effectM
}
