
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
##  Code tested on 2025-11-26


########################################################
## Functions to calculate p-values based on the 
##   asymptotic distribution: the tail probabilities
##   from the Square Ingegral of a Brownian Bridge
##
## We use the results of
## - Leonid Tolmatz. "On the Distribution of the Square 
##   Integral of the Brownian Bridge." Ann. Probab. 30 (1)
##   253 - 269, January 2002. https://doi.org/10.1214/aop/1020107767

integralSquaredBrownianBridgePvale <- function(x) {
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
  p_value <- integralSquaredBrownianBridgePvale(stat)
  c("SCUSUM Stat"=stat, "chpt Location"=k_change, "p-value"=p_value)
}


scusumTestProportion <- function(x) {
  ## x is a matrix, each row is #success, #trials
  n <- dim(x)[1]
  x <- as.matrix(x)
  
  no_change_phat <- sum(x[,1])/sum(x[,2])
  ## Wald estimate of standard error
  se_est <- sqrt(no_change_phat*(1-no_change_phat)/sum(x[,2]))
  
  ## Alternative approach is to take stdev of all sample p-hats
  ##   but Wald works better in simulations
  # se_est <- sd(x[,1]/x[,2])/sqrt(n)
  
  success_cusum <- cumsum(x[,1])
  trial_cusum <- cumsum(x[,2])
  
  regimetwo_success <- sum(x[,1]) - success_cusum
  regimetwo_total <- sum(x[,2]) - trial_cusum
  
  regime1_phats <- success_cusum/trial_cusum
  regime2_phats <- regimetwo_success/regimetwo_total
  regime2_phats[n] <- 0 # turn that NA into a 0
  
  k <- 1:n
  
  cusums.vals <- (k/n)*(1 - (k/n) )*((regime1_phats - regime2_phats))
  cusums.vals <- abs(cusums.vals)/se_est
  k_change <- which.max(abs(cusums.vals))
  stat <- mean(cusums.vals^2, na.rm=TRUE)
  if(is.na(stat)) {print(x)}
  p_value <- integralSquaredBrownianBridgePvale(stat)
  c("SCUSUM Stat"=stat, "chpt Location"=k_change, "p-value"=p_value)
}


