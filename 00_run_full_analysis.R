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

shorty_day_cutoff <- 3
source("01_fetching_best_track_data.R")
source("02_processing_hursat_adt.R")
rm(list=ls())
quarto::quarto_render("03_ibtracs_v_hursat.qmd")
source("04_tropical_cyclone_analysis.R")
rm(list=ls())
quarto::quarto_render("05_report_tropical_cyclone_analysis.qmd")
quarto::quarto_render("06_quick_updated_results.qmd")
source("07_results_ibtracs_plots.R")
rm(list=ls())
source("08_results_hursat_plots.R")
rm(list=ls())

