
################################################
##
##  04_tropical_cyclone_analysis.R
##
##  This file performs the analysis for the
##    tropical cyclone record.
##  After we input the processed data from the files
##      01_fetching_best_track_data.R
##      02_processing_hursat_adt.R
##  we then perform the analysis as outlined in the report.
##    Since the methods are the same for all basins
##    and data sources (IBTrACS or HURSAT) we put all
##    the analysis in a single function. We call the function
##    for the various combinations and save the result
##    to bundle into the technical report.
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-02-17
##      Ran on 2026-05-18 after adding in Hurricanes

library(tidyverse)
set.seed(1)

## Load data
load("./data/ibtracs_nonShortiesForProportions.RData")
load("./data/hursat_nonShortiesForProportions.RData")



## Load cusum test
source("cusumBasedChangePointsTest.R")

## Add the "Global" basin to each dataset
ibtracs_all_basins <- ibtracs_non_shorty_for_props |>
  group_by(SEASON) |>
  summarize(Hurricanes = sum(Hurricanes),
            Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms)) |>
  mutate(BASIN = "Global") |>
  bind_rows(ibtracs_non_shorty_for_props)

## Now with the HURSAT data
hursat_all_basins <- hursat_non_shorty_for_props |>
  group_by(SEASON) |>
  summarize(Hurricanes = sum(Hurricanes),
            Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms)) |>
  mutate(BASIN = "Global") |>
  bind_rows(hursat_non_shorty_for_props)

num_years_ibtracs <- length(unique(ibtracs_non_shorty_for_props$SEASON))
num_years_hursat <- length(unique(hursat_non_shorty_for_props$SEASON) )

## Trimmed version of IBTraCs matching HURSAT years
ibtracs_trimmed_all_basins <- ibtracs_all_basins |>
  dplyr::filter(SEASON %in% unique(hursat_all_basins$SEASON) )



## Main analysis function
##
## Specify a Basin and dataset (IBTracs vs HURSAT vs Trimmed IBTracs)
## All statistical models are put together
##    and provided in a named list.
get_report_basin <- function(basin="Global", all_basins=ibtracs_all_basins) {
  ## Some Data Processing for plots and summary statistics
  ##   and used in our analysis
  ## Note: YearsSinceStart is just a number counting the number of years
  ##   since the beginning of our record (1980 or 1981, depending on source)
  storm_basin <- all_basins |>
    dplyr::filter(BASIN==basin) |>
    mutate(YearsSinceStart = SEASON - min(SEASON)) |>
    mutate(Prop_Major = Major_Storms/Total_Storms,
           Prop_Intense = Intense_Storms/Total_Storms)
  
  ## Put data in a tall (or long) format for plotting
  storm_basin_long <- storm_basin |>
    dplyr::select(SEASON, Total_Storms, Hurricanes, Major_Storms, Prop_Major, Intense_Storms, Prop_Intense, BASIN) |>
    pivot_longer(c(Total_Storms, Hurricanes, Major_Storms, Prop_Major, Intense_Storms, Prop_Intense), 
                 names_to="Measure", values_to="values") |>
    mutate(Measure = factor(Measure, levels=c("Total_Storms", "Hurricanes", "Major_Storms", "Prop_Major", "Intense_Storms", "Prop_Intense"),
                            labels=c("Total Non-Shorties",
                                     "Hurricanes",
                                     "Major Storms",
                                     "Proportion of Major Storms",
                                     "Intense Storms",
                                     "Proportion of Intense Storms")))
  
  ## Time Series Plots of the data
  plot_data <- ggplot(storm_basin_long, aes(x=SEASON, y=values) ) +
    geom_line(color="gray50") + 
    geom_point(color="gray40") +
    facet_wrap(~Measure, nrow=5, scales="free_y") +
    theme_bw() + 
    theme(axis.title = element_blank() )
  
  ## Summary statistics (Full data aggregate)
  summary_stat <- storm_basin |>
    summarise(`Total Storms (Avg)` = mean(Total_Storms),
              `Total Storms (stddev)` = sd(Total_Storms),
              `Hurricanes (Avg)` = mean(Hurricanes),
              `Hurricanes (stddev)` = sd(Hurricanes),
              `Major Storms (Avg)` = mean(Major_Storms),
              `Major Storms (stddev)` = sd(Major_Storms),
              `Prop Major Storms` = sum(Major_Storms)/sum(Total_Storms),
              `Intense Storms (Avg)` = mean(Intense_Storms),
              `Intense Storms (stddev)` = sd(Intense_Storms),
              `Prop Intense Storms` = sum(Intense_Storms)/sum(Total_Storms) )
  
  ## Poisson Regression on Total, Hurricanes, Major & Intense
  
  count_pois_fit <- glm(Total_Storms ~ YearsSinceStart, data=storm_basin,
                        family=poisson(link="log"))
  
  hur_pois_fit <- glm(Hurricanes ~ YearsSinceStart, data=storm_basin,
                      family=poisson(link="log") )
  
  major_pois_fit <- glm(Major_Storms ~ YearsSinceStart, data=storm_basin,
                        family=poisson(link="log") )
  
  intense_pois_fit <- glm(Intense_Storms ~ YearsSinceStart, data=storm_basin,
                          family=poisson(link="log") )
  
  ## Poisson SCUSUM AMOC Changepoint Tests
  ##  on totals, hurricanes, major, and intense
  
  pois_chpt_test <- scusumTestPoisson(storm_basin$Total_Storms)
  
  hur_chpt_test <- scusumTestPoisson(storm_basin$Hurricanes)
  
  major_pois_chpt_test <- scusumTestPoisson(storm_basin$Major_Storms)
  
  intense_pois_chpt_test <- scusumTestPoisson(storm_basin$Intense_Storms)
  
  ## Logistic Regression on Proportion of Hurricanes,
  ##    Major and Proportion of Intense
  
  prop_data <- storm_basin |>
    mutate(Non_Hurricane = Total_Storms - Hurricanes,
           Non_Major = Total_Storms - Major_Storms,
           Non_Intense = Total_Storms - Intense_Storms)
  
  prop_hurr_logit_fit <-  glm(cbind(Hurricanes, Non_Hurricane) ~ YearsSinceStart, 
                              data=prop_data, 
                              family=binomial)
  
  prop_major_logit_fit <- glm(cbind(Major_Storms, Non_Major) ~ YearsSinceStart, 
                              data=prop_data, 
                              family=binomial)
  
  prop_intense_logit_fit <- glm(cbind(Intense_Storms, Non_Intense) ~ YearsSinceStart, 
                                data=prop_data, 
                                family=binomial)
  
  ## SCUSUM Test for Proportions
  ##    On Hurricanes, Major & Intense
  
  prop_hurr_chpt_test <- scusumTestProportion(cbind(prop_data$Hurricanes, prop_data$Total_Storms) )
  
  prop_major_chpt_test <- scusumTestProportion(cbind(prop_data$Major_Storms, prop_data$Total_Storms) )
  
  prop_intense_chpt_test <- scusumTestProportion(cbind(prop_data$Intense_Storms, prop_data$Total_Storms) )
  
  ## Output all results in a useful list
  list(basin_name = basin,
       summary_stat = summary_stat,
       
       pois_fit = count_pois_fit,
       hur_fit = hur_pois_fit,
       major_fit = major_pois_fit,
       intense_fit = intense_pois_fit,
       
       hurr_logit_fit = prop_hurr_logit_fit,
       major_logit_fit = prop_major_logit_fit,
       intense_logit_fit = prop_intense_logit_fit,
       
       pois_chpt_test = pois_chpt_test,
       hur_pois_chpt_test = hur_chpt_test,
       major_pois_chpt_test = major_pois_chpt_test,
       intense_pois_chpt_test = intense_pois_chpt_test,
       
       prop_hur_chpt_test = prop_hurr_chpt_test,
       prop_major_chpt_test = prop_major_chpt_test,
       prop_intense_chpt_test = prop_intense_chpt_test,
       
       plot_data = plot_data,
       
       data=storm_basin)
}

## Now call that function on all basins & data sources, saving
##   the results for the report.

ibtracs_global_results <- get_report_basin(basin="Global", all_basins=ibtracs_all_basins)
hursat_global_results <- get_report_basin(basin="Global", all_basins=hursat_all_basins)
ibtracs_trimmed_global_results <- get_report_basin(basin="Global", all_basins=ibtracs_trimmed_all_basins)

ibtracs_na_results <- get_report_basin(basin="North Atlantic", all_basins=ibtracs_all_basins)
hursat_na_results <- get_report_basin(basin="North Atlantic", all_basins=hursat_all_basins)
ibtracs_trimmed_na_results <- get_report_basin(basin="North Atlantic", all_basins=ibtracs_trimmed_all_basins)

ibtracs_ep_results <- get_report_basin(basin="Eastern North Pacific", all_basins=ibtracs_all_basins)
hursat_ep_results <- get_report_basin(basin="Eastern North Pacific", all_basins=hursat_all_basins)
ibtracs_trimmed_ep_results <- get_report_basin(basin="Eastern North Pacific", all_basins=ibtracs_trimmed_all_basins)

ibtracs_wp_results <- get_report_basin(basin="Western North Pacific", all_basins=ibtracs_all_basins)
hursat_wp_results <- get_report_basin(basin="Western North Pacific", all_basins=hursat_all_basins)
ibtracs_trimmed_wp_results <- get_report_basin(basin="Western North Pacific", all_basins=ibtracs_trimmed_all_basins)

ibtracs_ni_results <- get_report_basin(basin="Northern Indian", all_basins=ibtracs_all_basins)
hursat_ni_results <- get_report_basin(basin="Northern Indian", all_basins=hursat_all_basins)
ibtracs_trimmed_ni_results <- get_report_basin(basin="Northern Indian", all_basins=ibtracs_trimmed_all_basins)

ibtracs_si_results <- get_report_basin(basin="Southern Indian", all_basins=ibtracs_all_basins)
hursat_si_results <- get_report_basin(basin="Southern Indian", all_basins=hursat_all_basins)
ibtracs_trimmed_si_results <- get_report_basin(basin="Southern Indian", all_basins=ibtracs_trimmed_all_basins)

ibtracs_sp_results <- get_report_basin(basin="Southern Pacific", all_basins=ibtracs_all_basins)
hursat_sp_results <- get_report_basin(basin="Southern Pacific", all_basins=hursat_all_basins)
ibtracs_trimmed_sp_results <- get_report_basin(basin="Southern Pacific", all_basins=ibtracs_trimmed_all_basins)



save(ibtracs_global_results, ibtracs_na_results, ibtracs_ep_results,
     ibtracs_wp_results, ibtracs_ni_results, ibtracs_si_results, ibtracs_sp_results,
     num_years_ibtracs,
     file="./data/results_ibtracs.RData")
save(hursat_global_results, hursat_na_results, hursat_ep_results,
     hursat_wp_results, hursat_ni_results, hursat_si_results, hursat_sp_results,
     num_years_hursat,
     file="./data/results_hursat.RData")
save(ibtracs_trimmed_global_results, ibtracs_trimmed_na_results, ibtracs_trimmed_ep_results,
     ibtracs_trimmed_wp_results, ibtracs_trimmed_ni_results, ibtracs_trimmed_si_results, ibtracs_trimmed_sp_results,
     num_years_hursat,
     file="./data/results_trimmed.RData")



##################################################################
## Ultimately, we take all the results from above and
##   effectively summarize as one table in the manuscript
## This function helps build that table and streamlines
##   the code that builds the graphical results
##   where we plot the significant changepoint results.
## In the end, this code outputs in a tall format
##   - The Basin of interest (Global, NA, EP, WP, NI, SI, SP)
##   - The Method (Changepoint or Regression)
##   - The Series (TC, Hurricanes, Major, Prop_Major, Intense, Prop_Intense)
##   - The "Start" (Regime 1 for changepoint or predicted value at time 1 for regression)
##   - The "End" (Regime 2 for changepoint or predicted value at time N for regression)
##   - Shift - the observed shift for changepoint or the regression coefficient
##   - Relative_Shift - Percentage change for changepoint, to (End-Start)/N*10 -- linear slope at decadal level
##   - Time - the changepoint location (NA for Regression)
##   - p.value - the p-value for the relavent test (SCUSUM or LRT)

extract_test_results <- function(basin_result) {
  num_years <- length(unique(basin_result$data$SEASON))
  first_year <- min(unique(basin_result$data$SEASON))
  
  ## Build a data.frame with the Change points Results
  ##   Time Series, Change point location, p-value, Mean of Regime 1, Mean of Regime 2
  df1 <- data.frame(
    Response = c("Cyclone Freq", "Hurricanes", "Major Freq", "Prop Major", "Intense Freq", "Prop Intense"),
    `Chpt Location` = c(basin_result$pois_chpt_test[2]+first_year,
                        basin_result$hur_pois_chpt_test[2]+first_year,
                        basin_result$major_pois_chpt_test[2]+first_year,
                        basin_result$prop_major_chpt_test[2]+first_year,
                        basin_result$intense_pois_chpt_test[2]+first_year,
                        basin_result$prop_intense_chpt_test[2]+first_year),
    `p-value` =c(basin_result$pois_chpt_test[3],
                 basin_result$hur_pois_chpt_test[3],
                 basin_result$major_pois_chpt_test[3],
                 basin_result$prop_major_chpt_test[3],
                 basin_result$intense_pois_chpt_test[3],
                 basin_result$prop_intense_chpt_test[3] ),
    `Regime1` = c(mean(basin_result$data$Total_Storms[1:basin_result$pois_chpt_test[2]]),
                  mean(basin_result$data$Hurricanes[1:basin_result$hur_pois_chpt_test[2]]),
                  mean(basin_result$data$Major_Storms[1:basin_result$major_pois_chpt_test[2]]),
                  sum(basin_result$data$Major_Storms[1:basin_result$prop_major_chpt_test[2]])/sum(basin_result$data$Total_Storms[1:basin_result$prop_major_chpt_test[2]])*100,
                  mean(basin_result$data$Intense_Storms[1:basin_result$intense_pois_chpt_test[2]]),
                  sum(basin_result$data$Intense_Storms[1:basin_result$prop_intense_chpt_test[2]])/sum(basin_result$data$Total_Storms[1:basin_result$prop_intense_chpt_test[2]])*100 ),
    `Regime2` = c(mean(basin_result$data$Total_Storms[(basin_result$pois_chpt_test[2]+1):num_years]),
                  mean(basin_result$data$Hurricanes[(basin_result$hur_pois_chpt_test[2]+1):num_years]),
                  mean(basin_result$data$Major_Storms[(basin_result$major_pois_chpt_test[2]+1):num_years]),
                  sum(basin_result$data$Major_Storms[(basin_result$prop_major_chpt_test[2]+1):num_years])/sum(basin_result$data$Total_Storms[(basin_result$prop_major_chpt_test[2]+1):num_years])*100,
                  mean(basin_result$data$Intense_Storms[(basin_result$intense_pois_chpt_test[2]+1):num_years]),
                  sum(basin_result$data$Intense_Storms[(basin_result$prop_intense_chpt_test[2]+1):num_years])/sum(basin_result$data$Total_Storms[(basin_result$prop_intense_chpt_test[2]+1):num_years])*100)
  )
  
  df1 <- df1 |>
    mutate(Change = Regime2-Regime1,
           Perc_Change = Change/Regime1*100,
           Basin = basin_result$basin_name,
           Method = "Changepoint",
           Mean = unlist(basin_result$summary_stat[1,c(1,3,5,7,8,10)]) ) |>
    dplyr::select(Basin, Method, Series = Response, Mean, Start=Regime1, End=Regime2, Shift=Change, Relative_Shift=Perc_Change, Time=Chpt.Location, p.value)
  
  df2 <- data.frame(Basin = basin_result$basin_name,
             Method = rep("Regression", 6),
             Series = c("Cyclone Freq", "Hurricanes", "Major Freq", "Prop Major", "Intense Freq", "Prop Intense"),
             Shift = c(coef(basin_result$pois_fit)[2], 
                       coef(basin_result$hur_fit)[2],
                       coef(basin_result$major_fit)[2],
                       coef(basin_result$major_logit_fit)[2],
                       coef(basin_result$intense_fit)[2],
                       coef(basin_result$intense_logit_fit)[2]),
             `Start` = c(predict(basin_result$pois_fit, type="response")[1], 
                         predict(basin_result$hur_fit, type="response")[1],
                         predict(basin_result$major_fit, type="response")[1], 
                         predict(basin_result$major_logit_fit, type="response")[1]*100, 
                         predict(basin_result$intense_fit, type="response")[1],
                         predict(basin_result$intense_logit_fit, type="response")[1]*100 ),
             `End` = c(predict(basin_result$pois_fit, type="response")[num_years],
                       predict(basin_result$hur_fit, type="response")[num_years],
                       predict(basin_result$major_fit, type="response")[num_years],
                       predict(basin_result$major_logit_fit, type="response")[num_years]*100,
                       predict(basin_result$intense_fit, type="response")[num_years],
                       predict(basin_result$intense_logit_fit, type="response")[num_years]*100),
             Time = rep(NA, 6),
             p.value = c(drop1(basin_result$pois_fit, test="LRT")[2,5], drop1(basin_result$hur_fit, test="LRT")[2,5], 
                         drop1(basin_result$major_fit, test="LRT")[2,5], drop1(basin_result$major_logit_fit, test="LRT")[2,5], 
                         drop1(basin_result$intense_fit, test="LRT")[2,5], drop1(basin_result$intense_logit_fit, test="LRT")[2,5]) ) |>
    mutate(Magnitude =  (End- Start)/num_years_ibtracs*10,
           Mean = unlist(basin_result$summary_stat[1,c(1,3,5,7,8,10)]) ) |>
    dplyr::select(Basin, Method, Series, Mean, Start, End, Shift, Relative_Shift=Magnitude, Time, p.value)
  bind_rows(df1, df2)
}

ibtracs_key_results <- bind_rows(
  extract_test_results(ibtracs_global_results),
  extract_test_results(ibtracs_na_results),
  extract_test_results(ibtracs_ep_results),
  extract_test_results(ibtracs_wp_results),
  extract_test_results(ibtracs_ni_results),
  extract_test_results(ibtracs_si_results),
  extract_test_results(ibtracs_sp_results)
)


hursat_key_results <- bind_rows(
  extract_test_results(hursat_global_results),
  extract_test_results(hursat_na_results),
  extract_test_results(hursat_ep_results),
  extract_test_results(hursat_wp_results),
  extract_test_results(hursat_ni_results),
  extract_test_results(hursat_si_results),
  extract_test_results(hursat_sp_results)
)



save(ibtracs_key_results, hursat_key_results,
     file="./data/results_overview.RData")

