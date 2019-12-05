# endoSwitch

The R package estimates an endogenous switching regression model using full maximum likelihood estimation. 

The function can much replicate the results of the *movestay* command in STATA, though minor difference could occur due to differences in
the choice of optimization method.

To install the package, run the following functions in R:

install.packages("devtools") # Run this code if the devtools package is not installed.   
library(devtools)   
install_github("cbw1243/endoSwitch")  

The package is beta version. Further developments will include the calculation of treatment effects. 

Contact: Bowen Chen, PhD (bwchen@illinois.edu) 
