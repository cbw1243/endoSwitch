% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/endoSwitch2Stage.R
\name{endoSwitch2Stage}
\alias{endoSwitch2Stage}
\title{Endogenous switching regression}
\usage{
endoSwitch2Stage(data, OutcomeDep, SelectDep, OutcomeCov, SelectCov)
}
\arguments{
\item{data}{a data frame. Data for running the regression analysis.}

\item{OutcomeDep}{character. Dependent variable in the outcome equation.}

\item{SelectDep}{character. Dependent variable in the Selection model. The variable must be binary (0 or 1).}

\item{OutcomeCov}{character vector. Covariates in the outcome equation.}

\item{SelectCov}{character vector. Covariates in the selection equation.}
}
\value{
A list containing regression results.
}
\description{
This function estimates the endogenous switching regression model via two-stage estimations (Maddala, 1986)
}
\details{
The first stage uses a probit model to estimate the selection equation.
The second stage uses ordinary least squares including the inverse mills ratios computed from the first stage
estimation results to estimate the outcome equations.
}
\examples{

data(ImpactData)
OutcomeDep <- 'Output'
SelectDep <- 'CA'
OutcomeCov <- c('Age')
SelectCov <- c('Age', 'Perception')
Results <- endoSwitch2Stage(ImpactData, OutcomeDep, SelectDep, OutcomeCov, SelectCov)
# First stage regression results
summary(Results$FirstStageReg)
# Second stage regression results -- non-adopter
summary(Results$SecondStageReg.0)
# Second stage regression results -- adopter
summary(Results$SecondStageReg.1)
}
\references{
Maddala, G. S. “Limited-Dependent and Qualitative Variables in Econometrics.” Cambridge Books. Cambridge University Press, 1986.
}
