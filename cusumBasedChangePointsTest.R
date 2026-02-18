
################################################
##
##  cusumBasedChangePointsTest.R
##
##  This file contains the functions necessary
##    to perform the CUSUM-based change point 
##    test for Poisson & Binomial parameters
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-02-17


########################################################
## Functions to calculate p-values based on the 
##   asymptotic distribution: the tail probabilities
##   from the Square Integral of a Brownian Bridge
##
## We use the results of
## - Leonid Tolmatz. "On the Distribution of the Square 
##   Integral of the Brownian Bridge." Ann. Probab. 30 (1)
##   253 - 269, January 2002. https://doi.org/10.1214/aop/1020107767

integralSquaredBrownianBridgePvalue <- function(x) {
  if(x < 0.10) {
    1 - 2*sqrt(2/pi)*exp(-1/(8*x) )
  } else {
    2/pi/sqrt(pi*x)*exp(-pi*pi*x/2)
  }
}


##################################################################
##
## Here are the two call functions. In both, you simply pass in
##    a time series.
##
## in scusumTestPoisson() the average of the square of
##   the CUSUM statistic is calculated assuming a Poisson
##   marginal distribution. This suggest the variance be 
##   estimated with the mean, since in the Poisson mean=variance
##
## in scusumTestProportion() the time series needs to be a matrix
##   where each row is a pair (successes, trials). the average
##   of the square of the CUSUM statistic is calculated 
##   on the sample proportion -- here we use the Wald approximation
##   for the standard error of the parameter.

scusumTestPoisson <- function(x) {
  n <- length(x)
  k <- 1:n
  cusums.vals <- (cumsum(x) - k*mean(x))/(sqrt(mean(x))*sqrt(n) )
  k_change <- which.max(abs(cusums.vals))
  stat <- mean(cusums.vals^2)
  p_value <- integralSquaredBrownianBridgePvalue(stat)
  c("SCUSUM Stat"=stat, "chpt Location"=k_change, "p-value"=p_value)
}



scusumTestProportion <- function(x) {
  
  ## This function is part of the variance of the test statistic
  ##    For large values of lambda it is approximately 1/lambda
  H <- function(lambda) {
    lim <- 1:10000
    sum(exp(-lambda + lim*log(lambda) - log(lim) - lfactorial(lim) ) )
  }
  
  
  ## x is a matrix, each row is #success, #trials
  n <- dim(x)[1]
  x <- as.matrix(x)
  
  ## P_hat for each time point, and overall phat under H_0
  phat_t <- x[,1]/x[,2]
  phat_t[which(is.na(phat_t))] <- 0
  phat_null <- sum(x[,1])/sum(x[,2])
  
  ## Mean of total storms, for variance estimate
  lambda <- mean(x[,2])
  var_est <- phat_null^2*exp(-lambda)*(1-exp(-lambda)) + phat_null*(1- phat_null)*H(lambda)
  
  k <- 1:n
  cusums.vals <- (cumsum(phat_t) - (k/n)*sum(phat_t))/(sqrt(var_est)*sqrt(n) )
  
  k_change <- which.max(abs(cusums.vals))
  stat <- mean(cusums.vals^2)
  p_value <- integralSquaredBrownianBridgePvalue(stat)
  c("SCUSUM Stat"=stat, "chpt Location"=k_change, "p-value"=p_value)
}

