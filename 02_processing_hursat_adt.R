################################################
##
##  02_processing_hurdata_adt.R
##
##  This file processes the RData file that 
##    houses the raw HURDAT ADT data
##
##  The raw HURDATA data can be downloaded from the 
##    supplemental files NOAA NCEI. We use the 
##    script provided in the github repo:
##       https://github.com/tjfisher19/hursat_adt_in_R
##
##  This script opens that file and classifies
##    each storm. 
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-05-18


library(tidyverse)
load("./data/hursat_adt_raw.RData")

############################################################
############################################################
## We require the following to constants to be declared
##   If they are not declared elsewhere (for example, in
##       
##       
##   then we set these as the defaults.)

if(!exists("hursat_beginning_year") ) hursat_beginning_year <- 1981;
if(!exists("shorty_day_cutoff") ) shorty_day_cutoff <- 3;

#############################################
## There are 115 storms in HURSAT
##   with a differing SID than
##   in the IBTrACS data
## We use the file provided in the HURSAT
##   record to change the SIDS
##   So they match
## One storm in HURSAT, SID="2024291N14279"
##   maps to two storms in IBTrACS.  The storm
##   formed in the Gulf in October 2024 and
##   made landfall in Belize.  A day and a half
##   later it went into the Pacific and reformed
##   as a Eastern Pacific storm!  We have to handle
##   that particular mapping carefully.

## Get the Alternative SIDS
hursat_alt_id <- read_csv("./data/ADT-HURSAT_list_of_storm_id.csv",
                          col_types = "cccc") |>
  dplyr::select(SID = `ADT-HURSAT SID`, ALT_SID = `ALT SID IBTrACS` ) |>
  drop_na() |>
  filter(SID != "2024291N14279")

## A simple check and replace
hursat_adt_raw <- hursat_adt_raw |>
  mutate(SID = case_when(SID=="2024291N14279" & ymd(Date)<ymd("2024-10-20") ~ "2024293N17275",
                         SID=="2024291N14279" & ymd(Date)>=ymd("2024-10-20") ~ "2024296N13260",
                         .default = SID)) |>
  left_join(hursat_alt_id, by="SID") |>
  mutate(SID = ifelse(is.na(ALT_SID), SID, ALT_SID))


#####################################
## Assigning a basin
##
## Storms that form around Mexico could origin in either the 
##   Eastern Pacific or North Atlantic Basin 
##   (longitude between 260 and 280 degrees).
##   We cannot determine this with a simple check of Longitude
##   so we build a polygon that is the Pacific Ocean
##   just to the southwest of Mexico.
##
## The object below called "mex_ep_basin_points" defines
##   the region of the Pacific ocean southwest of Mexico
## It is a polygon that is fairly triangular. 
##  
## Essentially, if a storm formed in that region 
##   it is an Eastern Pacific storm
## If near Mexico, but not in this region, North Atlantic,
##   likely forming the Gulf.

library(sf)
mex_ep_basin_points <- rbind(c(260, 20),
                             c(260, 0),
                             c(280, 0),
                             c(280, 7),
                             c(260, 20))
mex_ep_poly <- st_polygon(list(mex_ep_basin_points))
mex_ep_sf <- st_sfc(mex_ep_poly, crs = 4326) # WGS84

###################################################
## Key Processing for basin determination here
##
## From the HURSAT documentation, the SID includes
##   most of the information needed.
##
##   Extract the hemisphere and origin Lat & Long
##   Use that information to find the Basin
##   of origin (we do this by comparing the 
##   Longitude) -- one caveat is the areas
##   Around Mexico -- I use some spatial methods
##   to check if the origin location is southwest
##   of Mexico (see the description above).
##


hursat_adt_raw <- hursat_adt_raw |>
  mutate(Hemi = str_sub(SID, 8, 8),
         Origin_Lat = as.numeric(str_sub(SID, 9, 10)),
         Origin_Lon = as.numeric(str_sub(SID, 11, 13)) )

df_sf <- st_as_sf(hursat_adt_raw[,c("Origin_Lon", "Origin_Lat")], coords = c(1,2), crs = 4326)
hursat_adt_raw$in_mex_ep_poly <- st_within(df_sf, mex_ep_sf, sparse=FALSE)[,1]


##################################################
## Now the key processing where each
##   storm is uniquely classified
##   based on maximum attained windpseed,
##   streak duration (at least # consecutive days
##   at tropical storm strength), basin and
##   season.

hursat_storms_fct <- hursat_adt_raw |> 
  ################
  ## Classify the basin based on the Hemisphere and Longitude
  ##   and for the spot around Mexico, if in that region
  ## Based on definitions in IBTrACS documentation
  ##   Page 26 of
  ##     https://www.ncei.noaa.gov/sites/default/files/2025-04/IBTrACS_version4r01_Technical_Details.pdf
  mutate(BASIN = case_when(Hemi=="N" & Origin_Lon <= 30 ~ "NA",    # West of Africa Coast
                           Hemi=="N" & Origin_Lon <= 100 ~ "NI",   # West of Myanmar & Thailand
                           Hemi=="N" & Origin_Lon <= 180 ~ "WP",   # West of Int. Date Line
                           Hemi=="N" & Origin_Lon <= 260 ~ "EP",   # West of California Coast
                           Hemi=="N" & Origin_Lon >= 280 ~ "NA",   # East of Florida
                           Hemi=="N" & in_mex_ep_poly ~ "EP",      # Mexico area, EP region
                           Hemi=="N" ~ "NA",                       # Mexico area, Gulf region (NA region)
                           Hemi=="S" & Origin_Lon <= 20 ~ "SA",    # West of Africa coast
                           Hemi=="S" & Origin_Lon <= 135 ~ "SI",   # West of Australia
                           Hemi=="S" & Origin_Lon <= 290 ~ "SP",   # West of South America
                           .default = "SA")                        # Default is South Atlantic (dropped)
  ) |>
  ##################################
  ## We need to create a date time variable
  ##    One of the times is reported as 25 hours...
  ## It is a clear error because the dates do not match up
  ##    so we fix it, hence the multple mutate() statements
  mutate(Date = ymd(Date),
         Hour = as.numeric(str_sub(Time, 1, 2)),
         Minute = as.numeric(str_sub(Time, 3, 4)),
         Second = as.numeric(str_sub(Time, 5, 6)) ) |>
  mutate(Date = as.Date(ifelse(Hour >= 24, Date+ days(1), Date) ),
         Hour = ifelse(Hour >= 24, Hour - 24, Hour) ) |>
  mutate(Year = year(Date),
         Day = day(Date)) |>
  mutate(DateTime = make_datetime(year=year(Date), month=month(Date), day=day(Date),
                                  hour = Hour, min=Minute, sec=Second),
         Year = year(Date)) |>
  ########################################
  ## Now more standard stuff, filter out
  ##    date before 1981 (need full satellite coverage)
  dplyr::filter(Year >= hursat_beginning_year) |>
  ########################################
  ##  Here, we determine the length
  ##   as per consecutive observations
  ##   where a storm is tropical storm strength
  drop_na(WindSpeed) |>
  mutate(TC = ifelse(WindSpeed > 34, "TC", "No")) |>
  group_by(SID, streak_id = consecutive_id(TC)) %>%
  mutate(streak_length = row_number()) |>
  dplyr::filter(TC=="TC") |>
  mutate(streak_duration = max(DateTime) - min(DateTime)) |>
  #########################################
  ## Now the key summarization. For each storm
  ##   we define the season, basin, its duration
  ##   whether it is a "shorty", tropical
  ##   cyclone classification based on Simpson-Saffir
  group_by(SID) %>% 
  summarize(N=n(),                               ## Number of observations
            SEASON=first(Year),                  ## Save the Season/year
            BASIN=first(BASIN),                  ## Save the BASIN
            Wind = max(WindSpeed, na.rm=TRUE),
            Duration = max(DateTime)-min(DateTime),
            Duration_Streak = max(streak_duration)) |>
  mutate(Duration_Class = ifelse(Duration < shorty_day_cutoff*24*60*60, "Shorty", "Non_Shorty"),
         Duration_Class_Streak = ifelse(Duration_Streak < shorty_day_cutoff*24*60*60, "Shorty", "Non_Shorty"),
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
##   the number of Hurricanes (Cat 1, 2, 3, 4, 5)
##   the number of Major storms (Cat 3, 4, 5)
##   the number of Intense storms (Cat 4, 5)
##   and total number (of non-shorties)
## From here, we can easily calculate proportions
##
## We also perform some extra processing where we input a 0
##   in any case when a season, basin, storm type 
##   combination did not occur (ie. Northern Indian)


hursat_non_shorty_for_props <- hursat_storms_fct |>
  dplyr::filter(Duration_Class_Streak=="Non_Shorty") |>
  group_by(SEASON, BASIN) |>
  summarize(Hurricanes = sum(StormCat !="Cat-0"),
            Major_Storms = sum(StormCat %in% c("Cat-3", "Cat-4", "Cat-5")),
            Intense_Storms = sum(StormCat %in% c("Cat-4", "Cat-5") ),
            Total_Storms = n() ) |>
  ungroup() |>
  complete(SEASON, BASIN, fill=list(Hurricanes=0, Major_Storms=0, Intense_Storms=0, Total_Storms=0)) %>% 
  dplyr::filter(BASIN != "Southern Atlantic")

save(hursat_non_shorty_for_props, 
     file="./data/hursat_nonShortiesForProportions.RData")


