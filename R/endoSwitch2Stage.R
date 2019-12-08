#' Endogenous switching regression
#'
#' This function estimates the endogenous switching regression model through two-stage regression.
#' The estimation procedures are described in Heckman et al. (2001). The \code{endoSwitch} function is
#' preferred than this function for estimation.
#'
#' The first stage uses a probit model to estimate the selection equation. The second stage uses OLS, plus the inverse mills ratios to estimate the main equations.
#'
#' @param RegData Data for doing the regression analysis
#' @param OutcomeDep Dependent variable in the main model
#' @param SelectDep Dependent variable in the selection model. This should be binary (0 or 1).
#' @param OutcomeCov Independent variables in the main model.
#' @param SelectCov Independent variable in the selection model
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
  MillsRatioAdt1 <- stats::dnorm(stats::predict(RegS1, RegData[AdtObs, ]))/stats::pnorm(stats::predict(RegS1, RegData[AdtObs, ]))
  MillsRatioAdt0 <- stats::dnorm(stats::predict(RegS1, RegData[-AdtObs, ]))/(1-stats::pnorm(stats::predict(RegS1, RegData[-AdtObs, ])))

  RegS2Formula_Adt1 <- stats::as.formula(paste0(OutcomeDep, '~', paste(OutcomeCov, collapse = '+'), '+', 'MillsRatioAdt1'))
  RegS2_Adt1 <- stats::lm(RegS2Formula_Adt1, data = RegData[AdtObs, ])

  RegS2Formula_Adt0 <- stats::as.formula(paste0(OutcomeDep, '~', paste(OutcomeCov, collapse = '+'), '+', 'MillsRatioAdt0'))
  RegS2_Adt0 <- stats::lm(RegS2Formula_Adt0, data = RegData[-AdtObs, ])

  Results <- list(FirstStageReg = RegS1,
                  SecondStageReg.0 = RegS2_Adt0,
                  SecondStageReg.1 = RegS2_Adt1)
  return(Results)
}
