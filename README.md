# Changes in the Tropical Cyclone Record

## Overview

This repository contains all the code and Quarto files (as well as the rendered version) for the statistical analysis accompanying the manuscript "Changepoint Analysis Confirms Decline in Observed Global Tropical Cyclone Frequency".

The Plots in the manuscript (including the supplemental plots for the HURSAT and Trimmed IBTRacs data) can be found in the [plots folder](plots/).

All the analysis can be found in the <a href="report_tropical_cyclone_analysis.html" target="_blank">Tropical Cyclone Analysis Report</a>.

## Repo Content

The repository contains the following:

Preliminary Code:

* `01_fetching_best_track_data.R` - R code that downloads the IBTraCS data and does some processing.
* `01_processing_hursat_adt.R` - R code that processes the downloaded the HURSAT files (See the code for more detail).
* `cusumBasedChangePointsTest.R` - R code that implements the SCUSUM based test for changes in the mean of a Poisson distributed time series and changes in the Proportion from a Binomial process.
* `tropical_cyclone_analysis.R` - R code that perform all analysis looking for changes in the Tropical Cyclone Record.

Statistical Analysis of Tropical Cyclone Record:

* `report_tropical_cyclone_analysis.qmd` - The main file that builds the "report" (a Quarto rendered website).
   + `report_analysis_code` - a folder containing the `qmd` files for each of the six basins and global record.
* `report_tropical_cyclone_anlysis.html` - The rendered HTML file of the report.
* `results_plots.R` - R code that will make the plots.
   + `plots` - a folder with the generated plots showing the change point segmentation.

Simulation Study Code:

* `simulationStudyCode.R` - R code that will perform a simulation study looking at the empirical size and power of the Generalized Linear Model (Poisson and Logistic regression) and the SCUSUM change point methods.
   + Results of the simulation study are in `./data/simulationStudyResults.RData`
* `report_simulation_study.qmd` - A fairly simple Quarto file that bundled all the simulation results into a short report. 
* `report_simulation_study.html` - The rendered HTML report.

## Installation and Requirements

All analysis were conducted and tested using R version 4.5.2 (2025-10-31 ucrt). 

The code was run on both a Windows v11 and Ubuntu 24.04.4 installations of R. 

The supplied code does rely on the `tidyverse` suite of packages in R for all data handling and plots. 

To replicate the data analysis results:

1. Download all .R files from this repository.
2. Execute all code in `01_fetching_best_track_data.R`
   + Some minor editing may be required for directory paths.
   + If using RStudio, you can click "Source" to run all code in the file.
   + Processed data is saved in the subdirectory `/data`.
3. Download the HURSAT data and execute all code in `01_processing_hursat_adt.R`.
   + Again, you may need to change the directory structure for file I/O.
   + Processed data is saved in the subdirectory `/data`.
4. Run all code in `tropical_cyclone_analysis.R`.
   + This code will load the changepoint methods in `cusumBasedChangePointsTest.R`.
   + The results of all these analyses will be saved in the `/data` subfolder.
5. Render the Quarto file `report_tropical_cyclone_analysis.qmd`.
   + This file will input all the Quarto files in the subfolder `/report_analysis_code`
   + A report in HTML will be generated summarizing the results.
6. Run the code in `results_ibtracs_plots.R`, `results_hursat_plots.R` and `results_trimmed_ibtracs_plots.R` to generate the plots provided in the manuscript and supplemental material.
   + All plots are saved as PNG files in the `/plots` subfolder.
   
**Optional Code:** If interest, you can also test the methods via simulation.

1. Run all the code in the `simulationStudyCode.R` file.
   + Running all this code will take upwards of 1 hour run time as it performs a large simulation looking to mimic the data we analyze.
2. Render the Quarto file `report_simulation_study.qmd` to generate a report summarizing the simulation results.
   + Generally you see all methods have good Type I error performance.
   + You can simulate power under different alternative hypotheses.




