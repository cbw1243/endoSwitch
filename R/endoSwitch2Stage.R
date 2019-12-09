#' Endogenous switching regression
#'
#' This function estimates the endogenous switching regression model through two-stage regression.
#' The estimation procedures are described in Heckman et al. (2001). The \code{endoSwitch} function is
#' preferred than this function for estimation.
#'
#' The first stage uses a probit model to estimate the selection equation. The second stage uses OLS, plus the inverse mills ratios to estimate the main equations.
#'
#' @param RegData an data frame. Data for running the regression analysis.
#' @param OutcomeDep character. Dependent variable in the outcome equation.
#' @param SelectDep character. Dependent variable in the Selection model. The iable must be binary (0 or 1).
#' @param OutcomeCov character vector. Independent iables in the outcome equation.
#' @param SelectCov character vector. Independent variables in the selection equation.
#'
#' @references Heckman, Tobias, and Vytlacil (2001). Four parameters of interest in the evaluation of
#' social programs. \emph{Southern Economic Journal}.

#' @return A list containing regression results.
#'
#' @export
#' @examples
#' data(ImpactData)
#' OutcomeDep <- 'Output'
#' SelectDep <- 'CA'
#' OutcomeCov <- c('Age')
#' SelectCov <- c('Age', 'Perception')
#' Results <- endoSwitch2Stage(ImpactData, OutcomeDep, SelectDep, OutcomeCov, SelectCov)
#' # First stage regression results
#' summary(Results$FirstStageReg)
#' # Second stage regression results -- non-adopter
#' summary(Results$SecondStageReg.0)
#' # Second stage regression results -- adopter
#' summary(Results$SecondStageReg.1)

endoSwitch2Stage <- function(RegData, OutcomeDep, SelectDep, OutcomeCov, SelectCov){
  RegData <- data.table::as.data.table(RegData)
  # Adopters
  AdtObs <- which(RegData[, SelectDep, with = F] == 1)

  # First stage estimation
  RegS1Formula <- stats::as.formula(paste0(SelectDep, '~', paste(SelectCov, collapse = '+')))

  RegS1 <- RegS1 <- tryCatch(
    { # Ignore warnings, mostly harmless
      suppressWarnings(stats::glm(RegS1Formula, data = RegData, family = stats::binomial(link = 'probit')))
    }, error = function(cond){
      message('Error in searching for starting values. Please provide start values manually.')
      message(paste0('The error message is: ', cond))
    }
  )

  # Second stage estimation
  Pred1 <- suppressWarnings(stats::predict(RegS1, RegData[AdtObs, ]))
  Pred0 <- suppressWarnings(stats::predict(RegS1, RegData[-AdtObs, ]))

  MillsRatioAdt1 <- stats::dnorm(Pred1)/stats::pnorm(Pred1)
  MillsRatioAdt0 <- stats::dnorm(Pred0)/(1-stats::pnorm(Pred0))

  RegS2Formula_Adt1 <- stats::as.formula(paste0(OutcomeDep, '~', paste(OutcomeCov, collapse = '+'), '+', 'MillsRatioAdt1'))
  RegS2_Adt1 <- stats::lm(RegS2Formula_Adt1, data = RegData[AdtObs, ])

  RegS2Formula_Adt0 <- stats::as.formula(paste0(OutcomeDep, '~', paste(OutcomeCov, collapse = '+'), '+', 'MillsRatioAdt0'))
  RegS2_Adt0 <- stats::lm(RegS2Formula_Adt0, data = RegData[-AdtObs, ])

  sigma1Est <- (sum(RegS2_Adt1$residuals^2) +
                    sum(MillsRatioAdt1*(MillsRatioAdt1 + Pred1))*stats::coef(RegS2_Adt1)['MillsRatioAdt1'])/length(MillsRatioAdt1)
  sigma1Est <- sigma1Est^.5

  sigma0Est <- (sum(RegS2_Adt0$residuals^2) +
                    sum(MillsRatioAdt0*(MillsRatioAdt0 + Pred0))*stats::coef(RegS2_Adt0)['MillsRatioAdt0'])/length(MillsRatioAdt0)
  sigma0Est <- sigma0Est^.5

  rho1Est <- stats::coef(RegS2_Adt1)['MillsRatioAdt1']/sigma1Est
  rho0Est <- -stats::coef(RegS2_Adt0)['MillsRatioAdt0']/sigma0Est

  endoSwitch2SResults <- list(FirstStageReg = RegS1, SecondStageReg.0 = RegS2_Adt0, SecondStageReg.1 = RegS2_Adt1,
                  distParEst = c(sigma0 = as.numeric(sigma0Est), sigma1 = as.numeric(sigma1Est),
                                 rho0 = as.numeric(rho0Est), rho1 = as.numeric(rho1Est)))
  return(endoSwitch2SResults)
}



