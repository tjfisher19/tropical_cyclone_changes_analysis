##############################################
## 
## 00_run_full_analysis.R
##
##  A convenience script that calls all
##    R code to run the analysis of the
##    Tropical Cyclone data.
##  Note, this scripe does not perform
##    the simulation study nor does it 
##    perform the shorty sensitivity
##    analysis. 
##
##  This is mainly here out of convenience
##    simply hit "source" and it should
##    run everything.
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-05-19


library(quarto)

## This variable defines the threshold for
##    so-called "shorties".  The code looks
##    at at least ## consecutive days
##    of 35+ knots intensity
shorty_day_cutoff <- 3

## Process the data
source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
rm(list=setdiff(ls(), "shorty_day_cutoff") )

## Report compaing data sources
quarto::quarto_render("03_ibtracs_v_hursat.qmd")

## Main analysis -- split across 3 files
##    The first file does all the modeling
##    The second generates a thorough report
##    The third provides summary tables
source("04_tropical_cyclone_analysis.R")
rm(list=setdiff(ls(), "shorty_day_cutoff") )
quarto::quarto_render("05_report_tropical_cyclone_analysis.qmd")
quarto::quarto_render("06_quick_updated_results.qmd")

## Make plots of changepoint results
source("07_results_ibtracs_plots.R")
rm(list=setdiff(ls(), "shorty_day_cutoff") )
source("08_results_hursat_plots.R")
rm(list=ls() )

