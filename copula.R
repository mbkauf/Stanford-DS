`" 
*-------------------------------------------------------------------------------
* Creating and Sampling from a Copula
* copula.R
* Date Created: 04/13/2022
* Last Updated: 
* Matthew Kaufmann, Stanford University
*-------------------------------------------------------------------------------
"`

## Load Packages
library(fitdistrplus)
library(MASS)
library(copula)
library(readxl)

## Parameters and Data to Build Copula
num  <- 1000
df   <- read_excel("rand_var.xlsx")

## Function
rand_copula <- function(seed=11111, n, data=df) {
  ## When examining the data, I determined that var1 should be fit to a gamma
  ## distribution and var2 should be fit to a beta distribution
  fit_var1 <- fitdist(data$var1, "gamma")
  fit_var2 <- fitdist(data$var2, "beta", method = "mme")
  cor_mat <- cor(data)  # Get correlation matrix between two variables
  
  set.seed(seed)
  ## First we need to set up a multivariate normal distribution with the 
  ## correlation between our 2 variables. Then we need to set up the margins to
  ## be the distributions for the underlying data
  myCop <- normalCopula(param=cor_mat[1,2], dim = 2, dispstr = "un") 
  myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta"),
                paramMargins=list(list(shape = fit_var1$estimate[1], 
                                       rate = fit_var1$estimate[2]),
                                  list(shape1 = fit_var2$estimate[1], 
                                       shape2 = fit_var2$estimate[2])))
  ## Now we sample from our copula and store sampled values as Z
  Z <<- rMvdc(n, myMvd)
  colnames(Z) <<- c("var1", "var2")
  Z[,2] <<- round(Z[,2],2)  # round var2 to 2 digits
}
rand_copula(n=num)
