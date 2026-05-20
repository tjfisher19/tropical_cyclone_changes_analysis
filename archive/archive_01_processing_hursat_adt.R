################################################
##
##  01_processing_hurdata_adt.R
##
##  This file processes the CSV files that 
##    house the HURDAT ADT data
##
##  The data can be downloaded from the 
##    supplemental files in the manuscript:
##
##  "Global increase in major tropical cyclone exceedance probability over the past four decades"
##  James P. Kossin, Kenneth R. Knapp, Timothy L. Olander, and Christopher S. Velden
##  Proc. Natl. Acad. Sci. U.S.A. 117 (22) 11975-11980
##  https://doi.org/10.1073/pnas.1920849117
##
##  We need 7 of the files (9 are supplied in the 
##    supplement) for our analysis
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2025-11-12

library(tidyverse)

#####################################################
## The data is spread across multiple files
##
## Each Row corresponds to a particular unique
##   tropical cyclone, and each column of the file
##   is the measurements associated with that
##   particular storms.
## So for example, the first row of data in each file
##   corresponds to the same tropical cyclone (storm 1978151N15260)
##   From the eastern Pacific Basin
##   The columns in the other file pairwise match
##   so the fifth 'wind' record matches the fifth year
##   and fifth month, fifth day and fifth time measurement.

storm_ids <- read_csv("./HURSAT_ADT/pnas.1920849117.sd07.csv")
basin <- read_csv("./HURSAT_ADT/pnas.1920849117.sd01.csv", na="")
wind <- read_csv("./HURSAT_ADT/pnas.1920849117.sd08.csv")
year_df <- read_csv("./HURSAT_ADT/pnas.1920849117.sd09.csv")
month_df <- read_csv("./HURSAT_ADT/pnas.1920849117.sd06.csv")
day_df <- read_csv("./HURSAT_ADT/pnas.1920849117.sd02.csv")
times_df <- read_csv("./HURSAT_ADT/pnas.1920849117.sd03.csv")

#############################################
## Link the 7 datasets into one by 
##   binding the relevant columns together
##   and pivoting from wide-to-tall format
## We also perform a little data processing

hursat_raw <- bind_cols(storm_ids, basin, year_df) |>
  pivot_longer(-c(StormID, Basin), values_to="Year") |>
  dplyr::select(-name) |>
  bind_cols(month_df |>
              pivot_longer(everything(), values_to="Month") |>
              dplyr::select(-name)) |>
  bind_cols(day_df |>
              pivot_longer(everything(), values_to="Day") |>
              dplyr::select(-name) ) |>
  bind_cols(times_df |>
              pivot_longer(everything(), values_to="Hour") |>
              dplyr::select(-name) ) |>
  bind_cols(wind |>
              pivot_longer(everything(), values_to="Wind") |>
              dplyr::select(-name) ) |>
  drop_na() |>
  mutate(Date = ymd(paste0(Year, "-", Month, "-", Day) ),
         DateTime = ymd_h(paste0(Year, "-", Month, "-", Day, "-", Hour)))

##################################################
## We perform a very similar processing as with
##   the IBTRacs data
## We remove all observations that are not at
##   Tropical Cyclone strength (at least 35 knots)
##   and start in 1981, the first complete year
##   of HURSAT data.
## Each unique storms is classified as a 
##   shorty or non-shorty based on duration
##   (48 or less hours, or more than 48 hours)
##    and Saffir-Simpson categorization (based
##    on maximum attained wind speed)

beginning_year <- 1981
current_year <- year(Sys.time())

hursat_storms_fct <- hursat_raw |>
  dplyr::filter(Wind > 34,
                Year >= beginning_year,
                Year < current_year) |>
  group_by(StormID) %>% 
  summarize(N=n(),                  ## Number of observations
            SEASON=first(Year),   ## Save the Season/year
            BASIN=first(Basin),     ## Save the BASIN
            Wind = max(Wind, na.rm=TRUE),
            Duration = max(DateTime)-min(DateTime)) |>
  mutate(Duration_Class = ifelse(Duration <= 48*60*60, "Shorty", "Non_Shorty"),
         BASIN = factor(BASIN,     ## Label the BASIN with something informative
                        levels=c("NA", "EP", "WP", "NI", "SI", "SP", "SA"),
                        labels=c("North Atlantic",
                                 "Eastern North Pacific",
                                 "Western North Pacific",
                                 "Northern Indian",
                                 "Southern Indian",
                                 "Southern Pacific",
                                 "Southern Atlantic")),
         StormCat = case_when(Wind >= 137 ~ "Cat-5",
                              Wind >= 113 ~ "Cat-4",
                              Wind >= 96  ~ "Cat-3",
                              Wind >= 83  ~ "Cat-2",
                              Wind >= 64  ~ "Cat-1",
                              Wind >= 34  ~ "Cat-0"), 
         StormCat = factor(StormCat, levels=c("Cat-0", "Cat-1", "Cat-2", "Cat-3", "Cat-4", "Cat-5") ) ) %>%
  dplyr::filter(!is.na(StormCat))

save(hursat_storms_fct, file = "./data/hursat_storm_classifications.RData")


############################################################
## Nnow build the dataset for our analysis.
##  We filter so only "non-shorties" are included
##  For each basin and season we count
##   the number of Major storms (Cat 3, 4, 5)
##   the number of Intense storms (Cat 4, 5)
##   and total number (of non-shorties)
## From here, we can easily calculate proportions
##
## We also perform some extra processing where we input a 0
##   in any case when a season, basin, storm type 
##   combination did not occur (ie. Northern Indian)


non_shorty_for_props <- hursat_storms_fct |>
  dplyr::filter(Duration_Class=="Non_Shorty") |>
  group_by(SEASON, BASIN) |>
  summarize(Major_Storms = sum(StormCat %in% c("Cat-3", "Cat-4", "Cat-5")),
            Intense_Storms = sum(StormCat %in% c("Cat-4", "Cat-5") ),
            Total_Storms = n() ) |>
  ungroup() |>
  complete(SEASON, BASIN, fill=list(Major_Storms=0, Intense_Storms=0, Total_Storms=0)) %>% 
  dplyr::filter(BASIN != "Southern Atlantic")

save(non_shorty_for_props, 
     file="./data/hursat_nonShortiesForProportions.RData")

