rm(list = ls())
data(ImpactData)
ManDepVar <- 'Output'
SelDepVar <- 'CA'
ManCovVar <- c('Age', 'Education')
SelCovVar <- c('Age', 'Education', 'Perception', 'Eastern', 'Southern')

Results <- endoSwitch(ImpactData, ManDepVar, SelDepVar, ManCovVar, SelCovVar)

summary(Results)
calcPar(Results)

remove.packages('endoSwitch')
devtools::install_github('cbw1243/endoSwitch')
library(endoSwitch)
?endoSwitch






