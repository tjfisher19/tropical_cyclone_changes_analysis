################################################
##
##  01_fetching_best_track_data.R
##
##  This file downloads the IBTracs data
##     and does all necessary processing
##     for the analysis.
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2025-11-12
##
##     Last ran on 2026-04-08
##         edited to use the USA windspeed records



## We use tidyverse to complete the processing

library(tidyverse)

################################
## Note, the first row of the csv file
##   is the variable names
## The second row is the units!
##   The data starts on line 3
## 
## So first read in the variable names only
##   then read in the data.
##
## Also note, the North Atlantic is coded as BASIN "NA"
##   which R reads as an NA - not available...
##

## Uncomment out the lines below to fetch the data
##    and save the result

# URL <- "https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.since1980.list.v04r01.csv"
# storm_col_names <- read_csv(URL,
#                             n_max = 0)
# 
# ibtracs_storms_raw <- read_csv(URL,
#                        skip = 2, col_names = names(storm_col_names), na="" )
# 
# save(ibtracs_storms_raw, file="./data/ibtracs_rawData.RData")

#######################################################
## Load Best Track Data
##   assuming you have previous fetched the raw data

load("./data/ibtracs_rawData.RData")

############################################################
############################################################
## We require the following to constants to be declared
##   If they are not declared elsewhere (for example, in
##       
##       
##   then we set these as the defaults.)

if(!exists("current_year") ) current_year <- year(Sys.time());
if(!exists("shorty_day_cutoff") ) shorty_day_cutoff <- 3;

##################################################################
## Each row of storms raw is a specific time & place observation
##   of a storm that is in the record. 
##
## For each storm, we do the following:
##  - Convert its SEASON into a numeric
##  - Calculate the maximum wind record across all the 
##    potential wind measurements for a given observation
##  - Remove all observations with wind speed less than 35 knots
##    since they are not tropical cyclones, and all storms
##    from the current (non-complete) calendar year.
##  - For each unique storm, we then find
##    + the number of observations in which it was 
##      a tropical cyclone
##    + Recording the SEASON
##    + Recording the BASIN
##    + Recording the name
##    + Recording the maximum attained wind speed (knots)
##    + Recording the duration of the storm
##  - We then classify each storm into a 
##    + Shorty - Lasting 48 hours or less
##    + Non_Shorty - Lasting more than 48 hours
##  - Add labels to the BASIN for plots
##  - Categorize each storm based on Saffir-Simpson scale
##
## When complete, each observation in storms_fct
##   will be a unique storm recording its season, basin
##   and classifications

ibtracs_storms_fct <- ibtracs_storms_raw %>%
  mutate(SEASON=year(ISO_TIME) ) |>
  mutate(Wind = USA_WIND ) |>
  drop_na(Wind) |>
  #mutate(Wind = pmax(WMO_WIND, USA_WIND, TOKYO_WIND, CMA_WIND, HKO_WIND, NEWDELHI_WIND, REUNION_WIND, BOM_WIND, NADI_WIND, WELLINGTON_WIND, DS824_WIND, TD9636_WIND, NEUMANN_WIND, MLC_WIND, na.rm=TRUE) ) |>
  # dplyr::filter(Wind > 34,
               # SEASON < current_year) |>
  dplyr::filter(SEASON < current_year,
                hour(ISO_TIME) %in% c(0,6,12,18),
                minute(ISO_TIME)==0)  |>
  mutate(TC = ifelse(Wind > 34, "TC", "No")) |>
  # mutate(new_streak = TC != lag(TC, default = first(TC)),
  #        streak_id = cumsum(new_streak) ) |>
  group_by(SID, streak_id = consecutive_id(TC)) %>%
  mutate(streak_length = row_number()) |>
  dplyr::filter(TC=="TC") |>
  group_by(SID) %>% 
  summarize(N=n(),                                 ## Number of observations
            SEASON=first(SEASON),                  ## Save the Season/year
            BASIN=first(BASIN),                    ## Save the BASIN of cycolgenesis
            NAME=first(NAME),                      ## Storm name (might be useful for annotation)
            Wind = max(Wind, na.rm=TRUE),          ## Maximum attained windspeed, for categorization
            Duration = max(ISO_TIME)-min(ISO_TIME), ## Duration of start and end of 35kn records
            Max_TC_Streak = max(streak_length)
            ) |>
  mutate(Duration_Class = ifelse(Duration < shorty_day_cutoff*24*60*60, "Shorty", "Non_Shorty"),
         Duration_Class_Streak = ifelse(Max_TC_Streak <= shorty_day_cutoff*4, "Shorty", "Non_Shorty"),
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

save(ibtracs_storms_fct, file = "./data/ibtracs_storm_classifications.RData")

############################################################
## We now build the dataset for our analysis.
##  We filter so only "non-shorties" are included
##  For each basin and season we count
##   the number of Hurricanes (Cat 1, 2, 3, 4, 5)
##   the number of Major storms (Cat 3, 4, 5)
##   the number of Intense storms (Cat 4, 5)
##   and total number (of non-shorties)
## From here, we can easily calculate proportions
##
## We also perform some extra processing where we input a 0
##   in any case when a season, basin, storm type 
##   combination did not occur (ie. Northern Indian)

ibtracs_non_shorty_for_props <- ibtracs_storms_fct |>
  dplyr::filter(Duration_Class_Streak=="Non_Shorty") |>
  group_by(SEASON, BASIN) |>
  summarize(Hurricanes = sum(StormCat !="Cat-0"),
            Major_Storms = sum(StormCat %in% c("Cat-3", "Cat-4", "Cat-5")),
            Intense_Storms = sum(StormCat %in% c("Cat-4", "Cat-5") ),
            Total_Storms = n()) |>
  ungroup() |>
  complete(SEASON, BASIN, fill=list(Hurricanes=0, Major_Storms=0, Intense_Storms=0, Total_Storms=0)) %>% 
  dplyr::filter(BASIN != "Southern Atlantic")

save(ibtracs_non_shorty_for_props, 
     file="./data/ibtracs_nonShortiesForProportions.RData")
