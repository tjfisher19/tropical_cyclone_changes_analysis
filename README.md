# Changes in the Tropical Cyclone Record

This repository contains all the code and Quarto files (as well as the rendered version) for the statistical analysis accompanying the manuscript "Changepoint Analysis Confirms Decline in Observed Global Tropical Cyclone Frequency".

The Plots in the manuscript (including the supplemental plots for the HURSAT and Trimmed IBTRacs data) can be found in the [plots folder](plots/).

All the analysis can be found in the [Tropical Cyclone Analysis Report](report_tropical_cyclone_analysis.html).

The repository contains the following:

Preliminary Code:

* `01_fetching_best_track_data.R` - R code that downloads the IBTracs data and does some processing.
* `01_processing_hursat_adt.R` - R code that processes the downloaded the HURSAT files (See the code for more detail).
* `cusumBasedChangePointsTest.R` - R code that implements the SCUSUM based test for changes in the mean of a Poisson distributed time series and changes in the Proportion from a Binomial process.

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

