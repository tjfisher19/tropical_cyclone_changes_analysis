################################################
##
##  20_study_shorty_sensitivity.R
##
##  Ultimately our analysis removes so-called
##    "shorties" from the record. These are
##    short duration storms that do not meet
##    the threshold to be a Tropical Cyclone
##
##  This code defines a "Shorty" based on
##    different duration (1 days vs 2.5 days)
##    and computes summary statistics and
##    looks at how it changes the statistical
##    analysis results.
##
##  After running this code, use
##    21_shorty_sensitivity_results.qmd
##    to generate a short report.
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-08-19

######################################
## Effectively, no shorties

shorty_day_cutoff <- 0

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key000 <- hursat_key_results
ibtracs_key000 <- ibtracs_key_results

hursat_global000 <- hursat_global_results
ibtracs_global000 <- ibtracs_global_results


#######################################
## At least 1 consecutive day of
##   35+ knots of wind speed

shorty_day_cutoff <- 1

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key010 <- hursat_key_results
ibtracs_key010 <- ibtracs_key_results

hursat_global010 <- hursat_global_results
ibtracs_global010 <- ibtracs_global_results

#######################################
## At least 1.5 consecutive days of
##   35+ knots of wind speed

shorty_day_cutoff <- 1.5

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key015 <- hursat_key_results
ibtracs_key015 <- ibtracs_key_results

hursat_global015 <- hursat_global_results
ibtracs_global015 <- ibtracs_global_results

#######################################
## At least 2 consecutive day2 of
##   35+ knots of wind speed

shorty_day_cutoff <- 2

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key020 <- hursat_key_results
ibtracs_key020 <- ibtracs_key_results

hursat_global020 <- hursat_global_results
ibtracs_global020 <- ibtracs_global_results


#######################################
## At least 2.5 consecutive days of
##   35+ knots of wind speed

shorty_day_cutoff <- 2.5

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key025 <- hursat_key_results
ibtracs_key025 <- ibtracs_key_results

hursat_global025 <- hursat_global_results
ibtracs_global025 <- ibtracs_global_results

#######################################
## At least 3 consecutive days of
##   35+ knots of wind speed

shorty_day_cutoff <- 3

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key030 <- hursat_key_results
ibtracs_key030 <- ibtracs_key_results

hursat_global030 <- hursat_global_results
ibtracs_global030 <- ibtracs_global_results


#######################################
## At least 3.5 consecutive days of
##   35+ knots of wind speed

shorty_day_cutoff <- 3.5

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key035 <- hursat_key_results
ibtracs_key035 <- ibtracs_key_results

hursat_global035 <- hursat_global_results
ibtracs_global035 <- ibtracs_global_results

#######################################
## At least 4 consecutive days of
##   35+ knots of wind speed

shorty_day_cutoff <- 4

source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
source("04_tropical_cyclone_analysis.R")

hursat_key040 <- hursat_key_results
ibtracs_key040 <- ibtracs_key_results

hursat_global040 <- hursat_global_results
ibtracs_global040 <- ibtracs_global_results



###############################
## Bundle all the resutls into
##   one data.frame

ibtracs_key_shorty <-
  bind_rows(
    ibtracs_key000 |> mutate(Days = 0),
    ibtracs_key010 |> mutate(Days = 1),
    ibtracs_key015 |> mutate(Days = 1.5),
    ibtracs_key020 |> mutate(Days = 2),
    ibtracs_key025 |> mutate(Days = 2.5),
    ibtracs_key030 |> mutate(Days = 3),
    ibtracs_key035 |> mutate(Days = 3.5),
    ibtracs_key040 |> mutate(Days = 4)
  )

hursat_key_shorty <-
  bind_rows(
    hursat_key000 |> mutate(Days = 0),
    hursat_key010 |> mutate(Days = 1),
    hursat_key015 |> mutate(Days = 1.5),
    hursat_key020 |> mutate(Days = 2),
    hursat_key025 |> mutate(Days = 2.5),
    hursat_key030 |> mutate(Days = 3),
    hursat_key035 |> mutate(Days = 3.5),
    hursat_key040 |> mutate(Days = 4)
  )


ibtracs_summaries <- 
  bind_rows(
    ibtracs_global000$summary_stat |> 
      mutate(Total = sum(ibtracs_global000$data$Total_Storms),
             Days = 0),
    ibtracs_global010$summary_stat |> 
      mutate(Total = sum(ibtracs_global010$data$Total_Storms),
             Days = 1),
    ibtracs_global015$summary_stat |> 
      mutate(Total = sum(ibtracs_global015$data$Total_Storms),
             Days = 1.5),
    ibtracs_global020$summary_stat |> 
      mutate(Total = sum(ibtracs_global020$data$Total_Storms),
             Days = 2),
    ibtracs_global025$summary_stat |> 
      mutate(Total = sum(ibtracs_global025$data$Total_Storms),
             Days = 2.5),
    ibtracs_global030$summary_stat |> 
      mutate(Total = sum(ibtracs_global030$data$Total_Storms),
             Days = 3),
    ibtracs_global035$summary_stat |> 
      mutate(Total = sum(ibtracs_global035$data$Total_Storms),
             Days = 3.5),
    ibtracs_global040$summary_stat |> 
      mutate(Total = sum(ibtracs_global040$data$Total_Storms),
             Days = 4)
  )


hursat_summaries <- 
  bind_rows(
    hursat_global000$summary_stat |> 
      mutate(Total = sum(hursat_global000$data$Total_Storms),
             Days = 0),
    hursat_global010$summary_stat |> 
      mutate(Total = sum(hursat_global010$data$Total_Storms),
             Days = 1),
    hursat_global015$summary_stat |> 
      mutate(Total = sum(hursat_global015$data$Total_Storms),
             Days = 1.5),
    hursat_global020$summary_stat |> 
      mutate(Total = sum(hursat_global020$data$Total_Storms),
             Days = 2),
    hursat_global025$summary_stat |> 
      mutate(Total = sum(hursat_global025$data$Total_Storms),
             Days = 2.5),
    hursat_global030$summary_stat |> 
      mutate(Total = sum(hursat_global030$data$Total_Storms),
             Days = 3),
    hursat_global035$summary_stat |> 
      mutate(Total = sum(hursat_global035$data$Total_Storms),
             Days = 3.5),
    hursat_global040$summary_stat |> 
      mutate(Total = sum(hursat_global040$data$Total_Storms),
             Days = 4)
  )

save(ibtracs_key_shorty, hursat_key_shorty,
     ibtracs_summaries, hursat_summaries,
     file="./data/shorty_sensitivity.RData")

