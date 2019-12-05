#' Endogenous switching regression
#'
#' This R function gives estimates of sigma and rho from an estimated model
#'
#'
#' @param Results An maxLik object estimated from endoSwitch
#'
#' @return A matrix that provides estimates, standard errors, and z values.
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
#'
calcPar <- function(Results){

  VarCov <- solve(-Results$hessian)
  coefEst <- coefficients(Results)

  SigmaNum <- grep('Sigma', names(coefEst))
  SigmaSD <- sqrt(diag(VarCov)[SigmaNum])

  RhoNum <- grep('Rho', names(coefEst))
  RhoSD <- sqrt(diag(VarCov)[RhoNum])

  outM <- cbind(Estimates = c(coefEst[SigmaNum], coefEst[RhoNum]),
                Std.error = c(SigmaSD, RhoSD))

  expFunc <- function(x) (exp(2*x) - 1)/(exp(2*x) + 1)

  parEst <- c(exp(outM[1, 1]), exp(outM[2, 1]), expFunc(outM[3, 1]), expFunc(outM[4, 1]))
  parSD <- c(msm::deltamethod (~ exp(x1), outM[1, 1], outM[1, 2]^2),
             msm::deltamethod (~ exp(x1), outM[2, 1], outM[2, 2]^2),
             msm::deltamethod (~ (exp(2*x1) - 1)/(exp(2*x1) + 1), outM[3, 1], outM[3, 2]^2),
             msm::deltamethod (~ (exp(2*x1) - 1)/(exp(2*x1) + 1), outM[4, 1], outM[4, 2]^2))

  parMatrix <- matrix(c(parEst, parSD, parEst/parSD), nrow = length(parEst))
  row.names(parMatrix) <- c('Main.0.Sigma', 'Main.1.Sigma', 'Main.0.Rho', 'Main.1.Rho')
  colnames(parMatrix) <- c('Estimate', 'Std. error', 'z')
  parMatrix

}
