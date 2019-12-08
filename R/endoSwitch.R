#' Endogenous switching regression
#'
#' This is the main interface for the endoSwitch package to estimate the endogenous switching regression
#' models.
#'
#' This function estimates the endogenous switching regression model using full maximum likelihood estimation
#' method. In this model, a selection equation sorts observation units over two different regimes (e.g., treated and
#' not-treated, or adopter or non-adopter), and two outcome equations that determine the
#' outcome. Estimation of the model relies on joint normality of the error terms in the equation system (the selection
#' equation plus two outcome equations). The model is estimated by maximizing the joint likelihood function that is provided in
#' Lokshin and Sajaia (2004).
#'
#' The \code{endoSwitch} uses the \code{maxLik} function in the maxLik package to do the optimization.
#'
#' Note that users can provide starting values for the optimization. Assume that you have M iables (including the constant) in
#' the selection equation, and N iables (including the constant) in the outcome equation.
#' Then you need (M + 2*N + 4) starting values. The first M values are for the variables in the
#' selection equation (last for the constant), then followed by N values for the outcome equation with SelectDep = 0,
#' and another N values for the outcome equation with SelectDep = 1.
#' The last four values are: sigma in outcome equation with SelectDep = 0, sigma in outcome equation with SelectDep = 1,
#' rho in outcome equation with SelectDep = 0, rho in outcome equation with SelectDep = 1. If starting values are not provided,
#' the function will automatically provide start values based on two-stage estimation results.
#'
#'
#' @param RegData an data frame. Data for running the regression analysis.
#' @param OutcomeDep character. Dependent variable in the outcome equation.
#' @param SelectDep character. Dependent variable in the Selection model. The iable must be binary (0 or 1).
#' @param OutcomeCov character vector. Independent iables in the outcome equation.
#' @param SelectCov character vector. Independent variables in the selection equation.
#' @param treatEffect TRUE/FALSE. Choose to show the average treatment effects or not.
#' @param trans TRUE/FALSE. Choose to show the estimates of original distributional parameters or not.
#' @param method character. Maximization method to be used. The default is "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno).
#' Other methods can also be used. See \code{\link{maxLik}}.
#' @param start optional numeric vector. Used as initial values of parameters for maximization purpose.
#' If NULL, the coefficient estimates from the two-stage estimation will be used.
#' @param verbose TRUE/FALSE. Choose to show the status of optimization or not.
#' @param ... Other parameters to be passed to the selected maximization routine. See \code{\link{maxLik}}.
#'
#' @return A list containing three elements. The first element is an object of class "maxLik". The estimates include parameters in the Selection equation,
#' parameters in the main equations, and the associated transformed distributional parameters. Note
#' that the distributional parameters have been transformed to facilitate the estimations, as recommended by Lokshin and Sajaia (2004).
#' As a default option, \code{endoSwitch} reports estimates of the distributional parameters by transforming them back to
#' origial forms using the Delta method. The second element contains the estimates of
#' distributional parameters.  The third element contains the treatment effects.
#'
#'
#' @references Lokshin and Sajaia (2004). Maximum-likelihood estimation of endogenous
#' switching regression models. \emph{The Stata Journal}.
#' HeckOutcome, Tobias, and Vytlacil (2001). Four parameters of interest in the evaluation of
#' social programs. \emph{Southern Economic Journal}.
#'
#' @export
#' @examples
#' data(ImpactData)
#' OutcomeDep <- 'Output'
#' SelectDep <- 'CA'
#' OutcomeCov <- c('Age')
#' SelectCov <- c('Age', 'Perception')
#' endoReg <- endoSwitch(ImpactData, OutcomeDep, SelectDep, OutcomeCov, SelectCov)
#' summary(endoReg$MLE.Results) # Sigma and Rho are transformed values.
#'
#' endoReg$distPar # Obtain the original distributional parameters
#'
#' endoReg$treatEffect # Obtain the treatment effect

endoSwitch <- function(RegData, OutcomeDep, SelectDep, OutcomeCov, SelectCov,
                       treatEffect = TRUE, trans = TRUE, method = 'BFGS', start = NULL,
                       verbose = FALSE, ...){
  if(sum(OutcomeDep %in% colnames(RegData)) != length(OutcomeDep))
    stop("OutcomeDep must be valid column names in the dataset.")

  if(sum(OutcomeCov %in% colnames(RegData)) != length(OutcomeCov))
    stop("OutcomeCov must be valid column names in the dataset.")

  if(sum(SelectDep %in% colnames(RegData)) != length(SelectDep))
    stop("OutcomeCov must be valid column names in the dataset.")

  if(sum(SelectCov %in% colnames(RegData)) != length(SelectCov))
    stop("OutcomeCov must be valid column names in the dataset.")

  if(nrow(RegData) < length(OutcomeCov) + length(SelectCov))
    stop('Too few observations.')

  if(identical(OutcomeCov, SelectCov))
    stop('OutcomeCov can not be identical to SelectCov.')

  RegData <- data.table::as.data.table(RegData)

  if(sum(apply(RegData[, c(OutcomeCov, OutcomeDep, SelectCov, SelectDep), with = FALSE],
               2, is.numeric)) < length(c(OutcomeCov, OutcomeDep, SelectCov, SelectDep)))
    stop('All Selected columns must be numeric.')

  SelectValue <- unlist(RegData[, SelectDep, with = FALSE])

  if(length(SelectValue[-c(which(SelectValue == 0), which(SelectValue == 1))]) >= 1)
    stop('The dependent iable (SelectDep) must be 0 or 1.')

  envir = new.env()
  environment(maxLogFunc) = envir

  ParNames <- c(paste0('Select.', SelectCov), 'Select.Const',
                paste0('Outcome.0.', OutcomeCov), 'Outcome.0.Const',
                paste0('Outcome.1.', OutcomeCov), 'Outcome.1.Const',
                'Outcome.0.SigmaX', 'Outcome.1.SigmaX',
                'Outcome.0.RhoX', 'Outcome.1.RhoX')

  assign('RegData', RegData, envir)
  assign('OutcomeDep', OutcomeDep, envir)
  assign('SelectDep', SelectDep, envir)
  assign('OutcomeCov', OutcomeCov, envir)
  assign('SelectCov', SelectCov, envir)

  assign('OutcomeParNum', length(OutcomeCov) + 1, envir)
  assign('SelectParNum', length(SelectCov) + 1, envir)
  assign('TotParNum', 2*(length(OutcomeCov) + 1) + length(SelectCov) + 1 + 4, envir)

  assign('verbose', verbose, envir)

  if(is.null(start)){
    twostage.results <- endoSwitch2Stage(RegData, OutcomeDep, SelectDep, OutcomeCov, SelectCov)
    start <- c(stats::coef(twostage.results$FirstStageReg)[-1], stats::coef(twostage.results$FirstStageReg)[1],
               stats::coef(twostage.results$SecondStageReg.0)[2:(length(OutcomeCov)+1)],
               stats::coef(twostage.results$SecondStageReg.0)[1],
               stats::coef(twostage.results$SecondStageReg.1)[2:(length(OutcomeCov)+1)],
               stats::coef(twostage.results$SecondStageReg.1)[1],
               1, 1, 0, 0)
  }

  if(length(start) != (length(SelectCov) + 1 + 2*(length(OutcomeCov) + 1) + 4))
    stop("wrong number of starting values")
  #if(isTRUE(verbose)) cat('Searching for optimal values.', '\r')

  mle.results <- maxLik::maxLik(maxLogFunc, start = start, method = method, ...)
  names(mle.results$estimate) <- names(mle.results$gradient) <- ParNames
  rownames(mle.results$hessian) <- colnames(mle.results$hessian) <- ParNames
  #if(isTRUE(verbose)) cat('Searching completed.', '\r')

  if(isTRUE(treatEffect)){
    treatEffectResult <- treatmentEffect(mle.results, RegData, OutcomeDep, SelectDep, OutcomeCov, SelectCov)
  }else{
    treatEffectResult <- 'Treatment effects are not calculated. Use treatEffect = TRUE to get them.'
  }

  if(isTRUE(trans)){
    distPar = calcPar(mle.results)
  }else{
    cat('Caution: Distributional parameters have been transformed. Use trans = TRUE to recover original estimates.')
    distPar <- NULL
  }

  out <- list(MLE.Results = mle.results, distPar = distPar, treatEffect = treatEffectResult)
  out
}
