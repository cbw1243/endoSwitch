#' Endogenous switching regression
#'
#' This function estimates the endogenous switching regression model through two-stage regression following Heckman et al. (2001).
#' The estimation procedures are described in
#'
#' The first stage uses a probit model to estimate the selection equation. The second stage uses OLS to estimate the main equations.
#'
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
#' # First stage regression results
#' summary(Results$FirstStageReg)
#' # Second stage regression results -- non-adopter
#' summary(Results$SecondStageReg_NoAdopt)
#' # Second stage regression results -- adopter
#' summary(Results$SecondStageReg_Adopt)

endoSwitch2Stage <- function(RegData, ManDepVar, SelDepVar, ManCovVar, SelCovVar){
  RegData <- data.table::as.data.table(RegData)
  # Adopters
  AdtObs <- which(RegData[, SelDepVar, with = F] == 1)

  # First stage estimation
  RegS1Formula <- as.formula(paste0(SelDepVar, '~', paste(SelCovVar, collapse = '+')))

  RegS1 <- stats::glm(RegS1Formula, data = RegData, family = binomial(link = 'probit'))

  # Second stage estimation
  MillsRatioAdt1 <- dnorm(predict(RegS1, ImpactData[AdtObs, ]))/pnorm(predict(RegS1, ImpactData[AdtObs, ]))
  MillsRatioAdt0 <- dnorm(predict(RegS1, ImpactData[-AdtObs, ]))/(1-pnorm(predict(RegS1, ImpactData[-AdtObs, ])))

  RegS2Formula_Adt1 <- as.formula(paste0(ManDepVar, '~', paste(ManCovVar, collapse = '+'), '+', 'MillsRatioAdt1'))
  RegS2_Adt1 <- lm(RegS2Formula_Adt1, data = RegData[AdtObs, ])

  RegS2Formula_Adt0 <- as.formula(paste0(ManDepVar, '~', paste(ManCovVar, collapse = '+'), '+', 'MillsRatioAdt0'))
  RegS2_Adt0 <- lm(RegS2Formula_Adt0, data = RegData[-AdtObs, ])

  Results <- list(FirstStageReg = RegS1,
                  SecondStageReg_NoAdopt = RegS2_Adt0,
                  SecondStageReg_Adopt =RegS2_Adt1)
  return(Results)
}
