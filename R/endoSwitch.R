#' Endogenous switching regression
#'
#' This is the main interface for the endoSwitch package, and the function performs Maximum Likelihood
#' estimation of an endogenous switching regression.
#'
#' This R function estimates parameters in an endogenous switching regression model by maximizing the
#' joint likelihood function, and the formula for calculating the likilihood function is from
#' Lokshin and Sajaia (2004).
#'
#' The function needs starting values to run. Assume that you have M variables (including the constant) in the selection model,
#' and N variables (including the constant) in the main model.
#' Then you need (M + 2*N + 4) starting values. The first M values are for the variables in the
#' selection model (last for the constant), the following N values for the main model with SelDepVar = 0,
#' and the following N values for the main model with SelDepVar = 1.
#' The last four values are: sigma in main model with SelDepVar = 0, sigma in main model with SelDepVar = 1,
#' rho in main model with SelDepVar = 0, rho in main model with SelDepVar = 1.
#'
#'
#' @param RegData Data for doing the regression analysis
#' @param ManDepVar Dependent variable in the main model
#' @param SelDepVar Dependent variable in the selection model. This should be binary (0 or 1).
#' @param ManCovVar Independent variables in the main model.
#' @param SelCovVar Independent variable in the selection model
#' @param method Maximization method.
#' @param start Numeric vector used as initial value of parameters for maximization purpose.
#' If NULL, the coefficient estimates from the two-stage estimation will be used.
#' @param verbose Print the log likelihood values during the optimization? TRUE or FALSE
#' @param ... Other parameters to be passed to the selected maximization routine.
#'
#' @return An object of class 'maxLik". The estimates include parameters in the selection model,
#' parameters in the main models in two different regimes, and the associated distributional parameters. Note
#' that the distributional parameters have been transformed to facilitate the estimations. Use calcPar function
#' to get estimates of original parameters.
#'
#' @export
#' @examples
#' data(ImpactData)
#' ManDepVar <- 'Output'
#' SelDepVar <- 'CA'
#' ManCovVar <- c('Age')
#' SelCovVar <- c('Age', 'Perception')
#' Results <- endoSwitch(ImpactData, ManDepVar, SelDepVar, ManCovVar, SelCovVar)
#' summary(Results)
#'
#' calcPar(Results)

endoSwitch <- function(RegData, ManDepVar, SelDepVar, ManCovVar, SelCovVar,
                       method = 'BFGS', start = NULL, verbose = FALSE, ...){
  RegData <- data.table::as.data.table(RegData)

  envir = new.env()
  environment(maxLogFunc) = envir

  ParNames <- c(paste0('Select.', SelCovVar), 'Select.Const',
                paste0('Main.0.', ManCovVar), 'Main.0.Const',
                paste0('Main.1.', ManCovVar), 'Main.1.Const',
                'Main.0.SigmaX', 'Main.1.SigmaX',
                'Main.0.RhoX', 'Main.1.RhoX')

  assign('RegData', RegData, envir)
  assign('ManDepVar', ManDepVar, envir)
  assign('SelDepVar', SelDepVar, envir)
  assign('ManCovVar', ManCovVar, envir)
  assign('SelCovVar', SelCovVar, envir)

  assign('ManParNum', length(ManCovVar) + 1, envir)
  assign('SelParNum', length(SelCovVar) + 1, envir)
  assign('TotParNum', 2*(length(ManCovVar) + 1) + length(SelCovVar) + 1 + 4, envir)

  assign('verbose', verbose, envir)

  if(is.null(start)){
    twostage.results <- endoSwitch2Stage(ImpactData, ManDepVar, SelDepVar, ManCovVar, SelCovVar)
    start <- c(coefficients(twostage.results$FirstStageReg)[-1], coefficients(twostage.results$FirstStageReg)[1],
              coefficients(twostage.results$SecondStageReg_NoAdopt)[2:(length(ManCovVar)+1)],
              coefficients(twostage.results$SecondStageReg_NoAdopt)[1],
              coefficients(twostage.results$SecondStageReg_Adopt)[2:(length(ManCovVar)+1)],
              coefficients(twostage.results$SecondStageReg_Adopt)[1],
              1, 1, 0, 0)
  }

  if(length(start) != (length(SelCovVar) + 1 + 2*(length(ManCovVar) + 1) + 4))
    stop("wrong number of starting values")

  mle.results <- maxLik::maxLik(maxLogFunc, start = start, method = method, ...)
  names(mle.results$estimate) <- names(mle.results$gradient) <- ParNames
  rownames(mle.results$hessian) <- colnames(mle.results$hessian) <- ParNames
  mle.results
}
