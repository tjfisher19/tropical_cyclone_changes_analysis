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
##  Code tested on 2026-01-15


source("cusumBasedChangePointsTest.R")

##########################
## Constants
runs=10000
n=length(1980:2025)

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
##  Total: 70.1087
##  Major: 24.69565
##    Prop Major: 0.3522481
##  Intense: 16.95652
##    Prop Intense: 0.2418605

global_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(70.1087, 70.1087))
global_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(24.69565, 24.69565))
global_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(70.1087, 70.1087), prop=c(0.3522481, 0.3522481))
global_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(16.95652, 16.95652))
global_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(70.1087, 70.1087), prop=c(0.2418605, 0.2418605))

global_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(global_total_null, global_major_null, global_major_prop_null, 
                                            global_intense_null, global_intense_prop_null))

global_null_reject_rates

#####################
## Alternative hypothesis

global_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 18, lambda=c(74.50000, 67.28571))
global_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 11, lambda=c(22.63636, 25.34286))
global_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 11, lambda=c(74.50000, 67.28571), prop=c(0.3055215, 0.3680498))
global_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 12, lambda=c(12.58333, 18.50000))
global_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(74.50000, 67.28571), prop=c(0.1712018, 0.2684592))

global_alt_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(global_total_alt, global_major_alt, global_major_prop_alt, 
                                            global_intense_alt, global_intense_prop_alt))

global_alt_reject_rates




##############################
## North Atlantic
#############
## Null hypothesis values (no changes)
##  Total: 10.45652
##  Major: 2.913044
##    Prop Major: 0.2785863
##  Intense: 1.847826
##    Prop Intense: 0.1767152

na_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(10.45652, 10.45652))
na_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.913044, 2.913044))
na_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(10.45652, 10.45652), prop=c(0.2785863, 0.2785863))
na_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(1.847826, 1.847826))
na_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(10.45652, 10.45652), prop=c(0.1767152, 0.1767152))

na_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(na_total_null, na_major_null, na_major_prop_null, 
                                            na_intense_null, na_intense_prop_null))

na_null_reject_rates

#####################
## Alternative hypothesis
##  

na_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 15, lambda=c(6.66667, 12.29032 ))
na_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 15, lambda=c(1.46667, 3.61290))
na_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 15, lambda=c(6.66667, 12.29032), prop=c(0.2200000, 0.2939633))
na_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 18, lambda=c(0.94444, 2.42857))
na_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 18, lambda=c(6.66667, 12.29032), prop=c(0.1268657, 0.1959654))

na_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(na_total_alt, na_major_alt, na_major_prop_alt, 
                                           na_intense_alt, na_intense_prop_alt))

na_alt_reject_rates





##############################
## Eastern Pacific
#############
## Null hypothesis values (no changes)
##  Total: 12.56522
##  Major: 4.673913
##    Prop Major: 0.376
##  Intense: 2.537
##    Prop Intense: 0.257

ep_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(12.56522, 12.56522))
ep_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(4.673913, 4.673913))
ep_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(12.56522, 12.56522), prop=c(0.3719723, 0.3719723))
ep_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(3.173913, 3.173913))
ep_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(12.56522, 12.56522), prop=c(0.2525952, 0.2525952))

ep_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(ep_total_null, ep_major_null, ep_major_prop_null, 
                                            ep_intense_null, ep_intense_prop_null))

ep_null_reject_rates

#####################
## Alternative hypothesis

ep_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 13, lambda=c(15.23077, 11.51515))
ep_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 19, lambda=c(5.42105, 4.14815))
ep_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 11, lambda=c(15.23077, 11.51515), prop=c(0.3291925, 0.3884892))
ep_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 34, lambda=c(2.64706, 4.66667))
ep_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(15.23077, 11.51515), prop=c(0.1666667, 0.2896040))

ep_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(ep_total_alt, ep_major_alt, ep_major_prop_alt, 
                                           ep_intense_alt, ep_intense_prop_alt))

ep_alt_reject_rates




##############################
## Western Pacific
#############
## Null hypothesis values (no changes)
##  Total: 21.93478
##  Major: 8.782609
##    Prop Major: 0.4003964
##  Intense: 6.804348
##    Prop Intense: 0.311

wp_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(21.93478, 21.93478))
wp_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(8.782609, 8.782609))
wp_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(21.93478, 21.93478), prop=c(0.4003964, 0.4003964))
wp_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(6.804348, 6.804348))
wp_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(21.93478, 21.93478), prop=c(0.3102081, 0.3102081))

wp_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(wp_total_null, wp_major_null, wp_major_prop_null, 
                                            wp_intense_null, wp_intense_prop_null))

wp_null_reject_rates

#####################
## Alternative hypothesis

wp_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 18, lambda=c(25.16667, 19.85714))
wp_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 37, lambda=c(9.16216, 7.22222))
wp_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 21, lambda=c(25.16667, 19.85714), prop=c(0.3598410, 0.4407115))
wp_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 7, lambda=c(4.71429, 7.17949))
wp_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 22, lambda=c(25.16667, 19.85714), prop=c(0.2604563, 0.3643892))

wp_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(wp_total_alt, wp_major_alt, wp_major_prop_alt, 
                                           wp_intense_alt, wp_intense_prop_alt))

wp_alt_reject_rates







##############################
## North Indian Basin
#############
## Null hypothesis values (no changes)
##  Total: 2.934783
##  Major: 0.8913043
##    Prop Major: 0.3037037
##  Intense: 0.6086957
##    Prop Intense: 0.2074074

ni_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.934783, 2.934783))
ni_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(0.8913043, 0.8913043))
ni_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(2.934783, 2.934783), prop=c(0.3037037, 0.3037037))
ni_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(0.6086957, 0.6086957))
ni_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(2.934783, 2.934783), prop=c(0.2074074, 0.2074074))

ni_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(ni_total_null, ni_major_null, ni_major_prop_null, 
                                            ni_intense_null, ni_intense_prop_null))

ni_null_reject_rates

#####################
## Alternative hypothesis

ni_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 12, lambda=c(2.25000, 3.17647))
ni_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 33, lambda=c(0.69697, 1.38462))
ni_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 26, lambda=c(2.25000, 3.17647), prop=c(0.2297297, 0.3934426))
ni_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 26, lambda=c(0.50000, 0.75000))
ni_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 26, lambda=c(2.25000, 3.17647), prop=c(0.1756757, 0.2459016))

ni_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(ni_total_alt, ni_major_alt, ni_major_prop_alt, 
                                           ni_intense_alt, ni_intense_prop_alt))

ni_alt_reject_rates




##############################
## Southern Indian
#############
## Null hypothesis values (no changes)
##  Total: 14.1087
##  Major: 4.956522
##    Prop Major: 0.3513097
##  Intense: 3.152174
##    Prop Intense: 0.2234206

si_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(14.1087, 14.1087))
si_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(4.956522, 4.956522))
si_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(14.1087, 14.1087), prop=c(0.3513097, 0.3513097))
si_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(3.152174, 3.152174))
si_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(14.1087, 14.1087), prop=c(0.2234206, 0.2234206))

si_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(si_total_null, si_major_null, si_major_prop_null, 
                                            si_intense_null, si_intense_prop_null))

si_null_reject_rates

#####################
## Alternative hypothesis

si_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 25, lambda=c(15.12000, 12.90476))
si_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 14, lambda=c(4.28571, 5.25000))
si_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 14, lambda=c(15.12000, 12.90476), prop=c(0.2912621, 0.3792325))
si_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 14, lambda=c(2.00000, 3.65625))
si_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 12, lambda=c(15.12000, 12.90476), prop=c(0.1195652, 0.2645161))

si_alt_reject_rates <- bind_cols(df_reject_rates,
                                     rbind(si_total_alt, si_major_alt, si_major_prop_alt, 
                                           si_intense_alt, si_intense_prop_alt))

si_alt_reject_rates





##############################
## South Pacific
#############
## Null hypothesis values (no changes)
##  Total: 8.108696
##  Major: 2.478261
##    Prop Major: 0.30563
##  Intense: 1.369565
##    Prop Intense: 0.1689008

sp_total_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(8.108696, 8.108696))
sp_major_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(2.478261, 2.478261))
sp_major_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(8.108696, 8.108696), prop=c(0.30563, 0.30563))
sp_intense_null <- sim_study_pois(runs=runs, n=n, chpt_loc = 20, lambda=c(1.369565, 1.369565))
sp_intense_prop_null <- sim_study_prop(runs=runs, n=n, chpt_loc = 20, lambda=c(8.108696, 8.108696), prop=c(0.1689008, 0.1689008))

sp_null_reject_rates <- bind_cols(df_reject_rates,
                                      rbind(sp_total_null, sp_major_null, sp_major_prop_null, 
                                            sp_intense_null, sp_intense_prop_null))

sp_null_reject_rates

#####################
## Alternative hypothesis

sp_total_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 21, lambda=c(10.09524, 6.44000))
sp_major_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 36, lambda=c(2.72222, 1.60000))
sp_major_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 23, lambda=c(10.09524, 6.44000), prop=c(0.2633929, 0.3691275))
sp_intense_alt <- sim_study_pois(runs=runs, n=n, chpt_loc = 23, lambda=c(1.08696, 1.65217))
sp_intense_prop_alt <- sim_study_prop(runs=runs, n=n, chpt_loc = 23, lambda=c(10.09524, 6.44000), prop=c(0.2550336, 0.1434264))

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



