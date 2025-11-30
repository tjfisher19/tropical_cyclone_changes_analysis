################################################
##
##  simulationStudyCode.R
##
##  This file contains the functions necessary
##    to perform simulation study required for
##    the simulation study report
##
##  Here, we simulate data to reasonably mimic
##    the observed tropical cyclone record
##    - When looking at counts, we simulated
##      Poisson data with a specified mean (lambda)
##    - For proportions, we simulate Poisson data
##      and a conditional Binomial, to then
##      calculate the proportion
##
##  The functions are set up to simulate either
##    under the null hypothesis (no changes)
##    or alternative hypothesis (a change occurs)
##    The method reports the rejection rate
##    of the GLM based method & SCUSUM based method
##    at an assortment of alpha-levels. It also
##    reports the "average" p-value and std.dev.
##
##  In the case of the "null" results, we want
##    rejection rates near the nominal level, and the
##    mean $p$-value should be about 0.50 with a standard
##    deviation of 0.28867 (sqrt(1/12)) according to standard 
##    mathematical statistical theory
##
##  In the case of the "alternative" hypothesis
##    the rejection rates should be higher than the
##    nominal rates, and the reported rejection rate
##    is an approximation of the "power" of the test
##    for the change point we found in the real data.
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2025-11-30


source("cusumBasedChangePointsTest.R")

##########################
## Constants
runs=10000
n=45
df_reject_rates <- data.frame(
  System = rep(c("Total Storms", "Major Storms", "Proportion Major", "Intense Storms", "Proportion Intense"),
               each=2),
  Method = rep(c("Regression", "CUSUM-Test"), 5)
)

###########################################
## Call functions to perform the study
sim_study_pois <- function(runs=1000, n=45, chpt_loc = 15, lambda=c(25,25) ) {
  
  lambda_vec <- c(rep(lambda[1], chpt_loc), rep(lambda[2], n-chpt_loc))
  N <- rep(n, runs)
  inside_fun <- function(n) {
    x <- rpois(n=n, lambda=lambda_vec)
    t <- (0:(n-1))
    count_pois_fit <- glm(x ~ t,
                          family=poisson(link="log"))
    c(drop1(count_pois_fit, test="LRT")$`Pr(>Chi)`[2], scusumTestPoisson(x)[3])
  }
  sample_pvals <- sapply(N, inside_fun)
  
  out <- rbind(c(`0.1% Rate` = mean(sample_pvals[1,] < 0.001)*100,
                 `1% Rate` = mean(sample_pvals[1,] < 0.01)*100,
                 `5% Rate` = mean(sample_pvals[1,] < 0.05)*100,
                 `10% Rate` = mean(sample_pvals[1,] < 0.10)*100,
                 `33% Rate` = mean(sample_pvals[1,] < 1/3)*100,
                 `Mean p-val` = mean(sample_pvals[1,]),
                 `SD p-val` = sd(sample_pvals[1,])),
               c(`0.1% Rate` = mean(sample_pvals[2,] < 0.001)*100,
                 `1% Rate` = mean(sample_pvals[2,] < 0.01)*100,
                 `5% Rate` = mean(sample_pvals[2,] < 0.05)*100,
                 `10% Rate` = mean(sample_pvals[2,] < 0.10)*100,
                 `33% Rate` = mean(sample_pvals[2,] < 1/3)*100,
                 `Mean p-val` = mean(sample_pvals[2,]),
                 `SD p-val` = sd(sample_pvals[2,]))
  )
  rownames(out) <- c("Poison Reg", "SCUSUM")
  out
}

sim_study_prop <- function(runs=1000, n=45, chpt_loc = 15, lambda=c(25,25), prop=c(0.3, 0.3) ) {
  
  lambda_vec <- c(rep(lambda[1], chpt_loc), rep(lambda[2], n-chpt_loc)) - 1
  prop_vec <- c(rep(prop[1], chpt_loc), rep(prop[2], n-chpt_loc))
  
  N <- rep(n, runs)
  inside_fun <- function(n) {
    x <- rpois(n=n, lambda=lambda_vec) + 1
    y <- rbinom(n=n, size=x, prob=prop_vec)
    t <- (0:(n-1))
    logit_fit <- glm(cbind(y, x-y) ~ t,
                     family=binomial(link="logit"))
    c(drop1(logit_fit, test="LRT")$`Pr(>Chi)`[2], scusumTestProportion(cbind(y,x))[3])
  }
  sample_pvals <- sapply(N, inside_fun)
  
  out <- rbind(c(`0.1% Rate` = mean(sample_pvals[1,] < 0.001)*100,
                 `1% Rate` = mean(sample_pvals[1,] < 0.01)*100,
                 `5% Rate` = mean(sample_pvals[1,] < 0.05)*100,
                 `10% Rate` = mean(sample_pvals[1,] < 0.10)*100,
                 `33% Rate` = mean(sample_pvals[1,] < 1/3)*100,
                 `Mean p-val` = mean(sample_pvals[1,]),
                 `SD p-val` = sd(sample_pvals[1,])),
               c(`0.1% Rate` = mean(sample_pvals[2,] < 0.001)*100,
                 `1% Rate` = mean(sample_pvals[2,] < 0.01)*100,
                 `5% Rate` = mean(sample_pvals[2,] < 0.05)*100,
                 `10% Rate` = mean(sample_pvals[2,] < 0.10)*100,
                 `33% Rate` = mean(sample_pvals[2,] < 1/3)*100,
                 `Mean p-val` = mean(sample_pvals[2,]),
                 `SD p-val` = sd(sample_pvals[2,]))
  )
  rownames(out) <- c("Logistic Reg", "SCUSUM")
  out
}


##############################
## Global "basin"
#############
## Null hypothesis values (no changes)
##  Total: 69.9
##  Major: 24.667
##    Prop Major: 0.353
##  Intense: 16.89
##    Prop Intense: 0.242

global_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(69.9, 69.9))
global_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(24.667, 24.667))
global_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(69.9, 69.9), prop=c(0.353, 0.353))
global_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(16.89, 16.89))
global_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(69.9, 69.9), prop=c(0.242, 0.242))

global_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(global_total_null, global_major_null, global_major_prop_null, 
                                            global_intense_null, global_intense_prop_null))

global_null_reject_rates

#####################
## Alternative hypothesis

global_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 18, lambda=c(74, 67.15))
global_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 12, lambda=c(22.58, 25.42))
global_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 11, lambda=c(73.82, 68.65), prop=c(0.304, 0.369))
global_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 12, lambda=c(12.42, 18.56))
global_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(72.83, 68.85), prop=c(0.17, 0.27))

global_alt_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(global_total_alt, global_major_alt, global_major_prop_alt, 
                                            global_intense_alt, global_intense_prop_alt))

global_alt_reject_rates




##############################
## North Atlantic
#############
## Null hypothesis values (no changes)
##  Total: 10.467
##  Major: 2.888
##    Prop Major: 0.276
##  Intense: 1.5
##    Prop Intense: 0.172

na_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(10.467, 10.467))
na_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.888, 2.888))
na_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(10.467, 10.467), prop=c(0.276, 0.276))
na_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(1.5, 1.5))
na_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(10.467, 10.467), prop=c(0.172, 0.172))

na_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(na_total_null, na_major_null, na_major_prop_null, 
                                            na_intense_null, na_intense_prop_null))

na_null_reject_rates

#####################
## Alternative hypothesis
##  

na_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 15, lambda=c(6.667, 12.367 ))
na_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 15, lambda=c(1.467, 3.6))
na_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 11, lambda=c(6.667, 12.367), prop=c(0.22, 0.29))
na_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 18, lambda=c(0.9444, 1.43))
na_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 18, lambda=c(6.667, 12.367), prop=c(0.127, 0.19))

na_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(na_total_alt, na_major_alt, na_major_prop_alt, 
                                           na_intense_alt, na_intense_prop_alt))

na_alt_reject_rates





##############################
## Eastern Pacific
#############
## Null hypothesis values (no changes)
##  Total: 12.467
##  Major: 4.688
##    Prop Major: 0.376
##  Intense: 2.537
##    Prop Intense: 0.257

ep_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(12.467, 12.467))
ep_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(4.688, 4.688))
ep_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(12.467, 12.467), prop=c(0.376, 0.376))
ep_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.537, 2.537))
ep_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(12.467, 12.467), prop=c(0.257, 0.257))

ep_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(ep_total_null, ep_major_null, ep_major_prop_null, 
                                            ep_intense_null, ep_intense_prop_null))

ep_null_reject_rates

#####################
## Alternative hypothesis

ep_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 13, lambda=c(15.23, 11.34))
ep_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 19, lambda=c(5.42, 4.15))
ep_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 34, lambda=c(15.23, 11.34), prop=c(0.358, 0.426))
ep_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 34, lambda=c(2.65, 4.9))
ep_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(15.23, 11.34), prop=c(0.218, 0.147))

ep_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(ep_total_alt, ep_major_alt, ep_major_prop_alt, 
                                           ep_intense_alt, ep_intense_prop_alt))

ep_alt_reject_rates




##############################
## Western Pacific
#############
## Null hypothesis values (no changes)
##  Total: 21.98
##  Major: 8.87
##    Prop Major: 0.403
##  Intense: 2.8
##    Prop Intense: 0.311

wp_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(21.98, 21.98))
wp_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(8.87, 8.87))
wp_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(21.98, 21.98), prop=c(0.403, 0.403))
wp_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.8, 2.8))
wp_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(21.98, 21.98), prop=c(0.311, 0.311))

wp_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(wp_total_null, wp_major_null, wp_major_prop_null, 
                                            wp_intense_null, wp_intense_prop_null))

wp_null_reject_rates

#####################
## Alternative hypothesis

wp_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 18, lambda=c(25.167, 19.85))
wp_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 37, lambda=c(9.162, 7.5))
wp_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 21, lambda=c(25.167, 19.85), prop=c(0.36, 0.45))
wp_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 7, lambda=c(4.71, 7.24))
wp_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 21, lambda=c(25.167, 19.85), prop=c(0.26, 0.37))

wp_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(wp_total_alt, wp_major_alt, wp_major_prop_alt, 
                                           wp_intense_alt, wp_intense_prop_alt))

wp_alt_reject_rates







##############################
## North Indian Basin
#############
## Null hypothesis values (no changes)
##  Total: 2.933
##  Major: 0.911
##    Prop Major: 0.310
##  Intense: 0.622
##    Prop Intense: 0.212

ni_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.933, 2.933))
ni_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(0.911, 0.911))
ni_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(2.933, 2.933), prop=c(0.310, 0.310))
ni_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(0.622, 0.622))
ni_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(2.933, 2.933), prop=c(0.212, 0.212))

ni_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(ni_total_null, ni_major_null, ni_major_prop_null, 
                                            ni_intense_null, ni_intense_prop_null))

ni_null_reject_rates

#####################
## Alternative hypothesis

ni_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 12, lambda=c(2.25, 3.18))
ni_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 33, lambda=c(0.70, 0.80))
ni_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 26, lambda=c(2.25, 3.18), prop=c(0.23, 0.41))
ni_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 26, lambda=c(0.5, 0.79))
ni_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 26, lambda=c(2.25, 3.18), prop=c(0.176, 0.26))

ni_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(ni_total_alt, ni_major_alt, ni_major_prop_alt, 
                                           ni_intense_alt, ni_intense_prop_alt))

ni_alt_reject_rates




##############################
## Southern Indian
#############
## Null hypothesis values (no changes)
##  Total: 13.89
##  Major: 4.8
##    Prop Major: 0.3456
##  Intense: 3.044
##    Prop Intense: 0.2192

si_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(13.89, 13.89))
si_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(4.8, 4.8))
si_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(13.89, 13.89), prop=c(0.3456, 0.3456))
si_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(3.044, 3.044))
si_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(13.89, 13.89), prop=c(0.2192, 0.2192))

si_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(si_total_null, si_major_null, si_major_prop_null, 
                                            si_intense_null, si_intense_prop_null))

si_null_reject_rates

#####################
## Alternative hypothesis

si_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 26, lambda=c(14.88, 12.53))
si_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 14, lambda=c(4.28, 5.03))
si_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(14.88, 12.53), prop=c(0.283, 0.371))
si_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 12, lambda=c(1.75, 3.52))
si_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(14.88, 12.53), prop=c(0.12, 0.26))

si_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(si_total_alt, si_major_alt, si_major_prop_alt, 
                                           si_intense_alt, si_intense_prop_alt))

si_alt_reject_rates





##############################
## South Pacific
#############
## Null hypothesis values (no changes)
##  Total: 8.18
##  Major: 2.5
##    Prop Major: 0.307
##  Intense: 1.38
##    Prop Intense: 0.168

sp_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(8.18, 8.18))
sp_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.5, 2.5))
sp_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(8.18, 8.18), prop=c(0.307, 0.307))
sp_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(1.38, 1.38))
sp_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(8.18, 8.18), prop=c(0.168, 0.168))

sp_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(sp_total_null, sp_major_null, sp_major_prop_null, 
                                            sp_intense_null, sp_intense_prop_null))

sp_null_reject_rates

#####################
## Alternative hypothesis

sp_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 21, lambda=c(10.09, 6.5))
sp_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 36, lambda=c(2.69, 1.77))
sp_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 23, lambda=c(10.09, 6.5), prop=c(0.261, 0.377))
sp_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 23, lambda=c(1.04, 1.72))
sp_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 23, lambda=c(10.09, 6.5), prop=c(0.108, 0.26))

sp_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(sp_total_alt, sp_major_alt, sp_major_prop_alt, 
                                           sp_intense_alt, sp_intense_prop_alt))

sp_alt_reject_rates



save(global_null_reject_rates, global_alt_reject_rates,
     na_null_reject_rates, na_alt_reject_rates,
     ep_null_reject_rates, ep_alt_reject_rates,
     wp_null_reject_rates, wp_alt_reject_rates,
     ni_null_reject_rates, ni_alt_reject_rates,
     si_null_reject_rates, si_alt_reject_rates,
     sp_null_reject_rates, sp_alt_reject_rates,
     file = "data/simulationStudyResults.RData")



