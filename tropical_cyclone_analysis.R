
################################################
##
##  tropical_cyclone_analysis.R
##
##  This file performs the analysis for the
##    tropical cyclone record.
##  After we input the processed data from the files
##      01_fetching_best_track_data.R
##      01_processing_hursat_adt.R
##  we then perform the analysis as outlined in the report.
##    Since the methods are the same for all basins
##    and data sources (IBTRaCS or HURSAT) we put all
##    the analysis in a single function. We call the function
##    for the various combinations and save the result
##    to bundle into the techincal report.
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-02-17

library(tidyverse)
set.seed(1)

## Load data
load("./data/ibtracs_nonShortiesForProportions.RData")
ibtracs_data <- non_shorty_for_props
load("./data/hursat_nonShortiesForProportions.RData")
hursat_data <- non_shorty_for_props

## Load cusum test
source("cusumBasedChangePointsTest.R")

## Add the "Global" basin to each dataset
ibtracs_all_basins <- ibtracs_data |>
  group_by(SEASON) |>
  summarize(Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms)) |>
  mutate(BASIN = "Global") |>
  bind_rows(ibtracs_data)

## Now with the HURSAT data
hursat_all_basins <- hursat_data |>
  group_by(SEASON) |>
  summarize(Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms)) |>
  mutate(BASIN = "Global") |>
  bind_rows(hursat_data)

num_years_ibtracs <- length(unique(ibtracs_data$SEASON))
num_years_hursat <- length(unique(hursat_data$SEASON) )

## Trimmed version of IBTraCs matching HURSAT years
ibtracs_trimmed_all_basins <- ibtracs_all_basins |>
  dplyr::filter(SEASON %in% unique(hursat_all_basins$SEASON) )



## Main analysis function
##
## Specify a Basin and dataset (IBTracs vs HURSAT vs Trimmed IBTracs)
## All analysis we need is reported.
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
    dplyr::select(SEASON, Total_Storms, Major_Storms, Prop_Major, Intense_Storms, Prop_Intense, BASIN) |>
    pivot_longer(c(Total_Storms, Major_Storms, Prop_Major, Intense_Storms, Prop_Intense), names_to="Measure", values_to="values") |>
    mutate(Measure = factor(Measure, levels=c("Total_Storms", "Major_Storms", "Prop_Major", "Intense_Storms", "Prop_Intense"),
                            labels=c("Total Non-Shorties",
                                     "Major Storms",
                                     "Proportion of Major Storms",
                                     "Intense Storms",
                                     "Proportion of Intense Storms")))
  
  ## Basic Plots of the data
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
              `Major Storms (Avg)` = mean(Major_Storms),
              `Major Storms (stddev)` = sd(Major_Storms),
              `Prop Major Storms` = sum(Major_Storms)/sum(Total_Storms),
              `Intense Storms (Avg)` = mean(Intense_Storms),
              `Intense Storms (stddev)` = sd(Intense_Storms),
              `Prop Intense Storms` = sum(Intense_Storms)/sum(Total_Storms) )
  
  ## Poisson Regression on Total, Major & Intense
  
  count_pois_fit <- glm(Total_Storms ~ YearsSinceStart, data=storm_basin,
                        family=poisson(link="log"))
  
  major_pois_fit <- glm(Major_Storms ~ YearsSinceStart, data=storm_basin,
                        family=poisson(link="log") )
  
  intense_pois_fit <- glm(Intense_Storms ~ YearsSinceStart, data=storm_basin,
                          family=poisson(link="log") )
  
  ## Poisson SCUSUM AMOC Changepoint Test
  ##  on totals, major, and intense
  
  pois_chpt_test <- scusumTestPoisson(storm_basin$Total_Storms)
  
  major_pois_chpt_test <- scusumTestPoisson(storm_basin$Major_Storms)
  
  intense_pois_chpt_test <- scusumTestPoisson(storm_basin$Intense_Storms)
  
  ## Logistic Regression on Proportion of Major
  ##    and Proportion of Intense
  
  prop_data <- storm_basin |>
    mutate(Non_Major = Total_Storms - Major_Storms,
           Non_Intense = Total_Storms - Intense_Storms)
  
  prop_major_logit_fit <- glm(cbind(Major_Storms, Non_Major) ~ YearsSinceStart, 
                              data=prop_data, 
                              family=binomial)
  
  prop_intense_logit_fit <- glm(cbind(Intense_Storms, Non_Intense) ~ YearsSinceStart, 
                                data=prop_data, 
                                family=binomial)
  
  ## SCUSUM Test for Proportions
  ##    On Major & Intense
  
  prop_major_chpt_test <- scusumTestProportion(cbind(prop_data$Major_Storms, prop_data$Total_Storms) )
  
  prop_intense_chpt_test <- scusumTestProportion(cbind(prop_data$Intense_Storms, prop_data$Total_Storms) )
  
  ## Output all results in a useful list
  list(summary_stat = summary_stat,
       
       pois_fit = count_pois_fit,
       major_fit = major_pois_fit,
       intense_fit = intense_pois_fit,
       
       major_logit_fit = prop_major_logit_fit,
       intense_logit_fit = prop_intense_logit_fit,
       
       pois_chpt_test = pois_chpt_test,
       major_pois_chpt_test = major_pois_chpt_test,
       intense_pois_chpt_test = intense_pois_chpt_test,
       
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



