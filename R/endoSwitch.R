#' Endogenous switching regression
#'
#' This is the main interface for the endoSwitch package, and the function performs full maximum likelihood
#' estimation of an endogenous switching regression.
#'
#' This function estimates the endogenous switching regression model using full maximum likelihood estimation
#' method. In this model, a selection equation sorts observation units over two different regimes (e.g., treated and
#' not-treated, adopter or non-adopter), and two main equations (also called outcome equations) determine the
#' outcome. The estimation of the model relies on joint normality of the error terms in the equation system (the selection
#' equation + 2 main equations). The model is estimated by maximizing the joint likelihood function, which is provided in
#' Lokshin and Sajaia (2004).
#'
#' The \code{endoSwitch} uses the \code{maxLik} function in the maxLik package to do the optimization.
#'
#' Note that users can provide starting values for the optimization. Assume that you have M variables (including the constant) in the selection equation,
#' and N variables (including the constant) in the main equation.
#' Then you need (M + 2*N + 4) starting values. The first M values are for the variables in the
#' selection equation (last for the constant), then followed by N values for the main model with SelDepVar = 0,
#' and another N values for the main model with SelDepVar = 1.
#' The last four values are: sigma in main equation with SelDepVar = 0, sigma in main equation with SelDepVar = 1,
#' rho in main equation with SelDepVar = 0, rho in main equation with SelDepVar = 1. If starting values are not provided,
#' the function will automatically search for reasonal values (usually works).
#'
#'
#' @param RegData an data frame. Data for running the regression analysis.
#' @param ManDepVar character. Dependent variable in the main model.
#' @param SelDepVar character. Dependent variable in the selection model. The variable must be binary (0 or 1).
#' @param ManCovVar character vector. Independent variables in the main model.
#' @param SelCovVar character vector. Independent variable in the selection model.
#' @param method character. Maximization method to be used. The default is "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno).
#' Other methods can also be used. See \code{\link{maxLik}}.
#' @param start optional numeric vector. Used as initial values of parameters for maximization purpose.
#' If NULL, the coefficient estimates from the two-stage estimation will be used.
#' @param verbose TRUE/FALSE. Showing the status of the optimization function.
#' @param ... Other parameters to be passed to the selected maximization routine.
#'
#' @return An object of class "maxLik". The estimates include parameters in the selection equation,
#' parameters in the main equations, and the associated distributional parameters. Note
#' that the distributional parameters have been transformed to facilitate the estimations. Use \code{calcPar} function
#' to get estimates of original parameters.
#'
#' @references Lokshin and Sajaia (2004). Maximum-likelihood estimation of endogenous
#' switching regression models. \emph{The Stata Journal}.
#' Heckman, Tobias, and Vytlacil (2001). Four parameters of interest in the evaluation of
#' social programs. \emph{Southern Economic Journal}.
#'
#' @export
#' @examples
#' data(ImpactData)
#' ManDepVar <- 'Output'
#' SelDepVar <- 'CA'
#' ManCovVar <- c('Age')
#' SelCovVar <- c('Age', 'Perception')
#' Results <- endoSwitch(ImpactData, ManDepVar, SelDepVar, ManCovVar, SelCovVar)
#' summary(Results) # Sigma and Rho are transformed values.
#'
#' calcPar(Results) # Obtain the original parameters

endoSwitch <- function(RegData, ManDepVar, SelDepVar, ManCovVar, SelCovVar,
                       method = 'BFGS', start = NULL, verbose = FALSE, ...){
  if(sum(ManDepVar %in% colnames(RegData)) != length(ManDepVar))
    stop("ManDepVar must be valid column names in the dataset.")

  if(sum(ManCovVar %in% colnames(RegData)) != length(ManCovVar))
    stop("ManCovVar must be valid column names in the dataset.")

  if(sum(SelDepVar %in% colnames(RegData)) != length(SelDepVar))
    stop("ManCovVar must be valid column names in the dataset.")

  if(sum(SelCovVar %in% colnames(RegData)) != length(SelCovVar))
    stop("ManCovVar must be valid column names in the dataset.")

  if(nrow(RegData) < length(ManCovVar) + length(SelCovVar))
    stop('Too few observations.')

  if(identical(ManCovVar, SelCovVar))
    stop('ManCovVar can not be identical to SelCovVar.')

  RegData <- data.table::as.data.table(RegData)

  if(sum(apply(RegData[, c(ManCovVar, ManDepVar, SelCovVar, SelDepVar), with = FALSE],
               2, is.numeric)) < length(c(ManCovVar, ManDepVar, SelCovVar, SelDepVar)))
    stop('All selected columns must be numeric.')

  SelValue <- unlist(RegData[, SelDepVar, with = FALSE])

  if(length(SelValue[-c(which(SelValue == 0), which(SelValue == 1))]) >= 1)
    stop('The dependent variable (SelDepVar) must be 0 or 1.')

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
    start <- c(stats::coef(twostage.results$FirstStageReg)[-1], stats::coef(twostage.results$FirstStageReg)[1],
               stats::coef(twostage.results$SecondStageReg_NoAdopt)[2:(length(ManCovVar)+1)],
               stats::coef(twostage.results$SecondStageReg_NoAdopt)[1],
               stats::coef(twostage.results$SecondStageReg_Adopt)[2:(length(ManCovVar)+1)],
               stats::coef(twostage.results$SecondStageReg_Adopt)[1],
               1, 1, 0, 0)
  }

  if(length(start) != (length(SelCovVar) + 1 + 2*(length(ManCovVar) + 1) + 4))
    stop("wrong number of starting values")
  if(isTRUE(verbose)) cat('Searching for optimal values.', '\r')

  mle.results <- maxLik::maxLik(maxLogFunc, start = start, method = method, ...)
  names(mle.results$estimate) <- names(mle.results$gradient) <- ParNames
  rownames(mle.results$hessian) <- colnames(mle.results$hessian) <- ParNames
  if(isTRUE(verbose)) cat('Searching completed.', '\r')
  mle.results
}
