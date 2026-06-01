################################################
##
##  results_ibtracs_plots.R
##
##  This file builds the plots of the change
##    point results. If the found p-value
##    is less than 1/3 (the "likely" finding)
##    we plot the suggested segmentation for
##    the Full IBTraCS data.
##
##  The first plot is all the findings at the
##    global level. Then the subsequent plots
##    are the results within each of the 6 basins
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-05-19


library(tidyverse)
library(patchwork)

rm(list=ls())
yearFirst <- 1980
yearLast <- 2025
sig_level <- 1/3


#######################################
## Load data & Change point findings

load("./data/ibtracs_nonShortiesForProportions.RData")
load("./data/results_overview.RData")

#########################################
## PLOT Parameters and
##   customization of axes and such
##   

if(!exists("shorty_day_cutoff") ) shorty_day_cutoff <- 3;

data_dots <- "gray50"
data_line <- "gray60"
seg_line <- "gray20"
caption_text <- "Source: International Best Track Archive for Climate Stewardship (IBTrACS)\nhttps://www.ncei.noaa.gov/products/international-best-track-archive"
subtitle_text <- paste("Segmentations based on significant SCUSUM Change Point Test\nStorms with less than", shorty_day_cutoff, "consecutive days at tropical cyclone strength excluded")

our_labels <- function(x) {
  case_when(x==yearFirst ~ paste0(yearFirst),
            x==2025 ~ "2025",
            x==2000 ~ "2000",
            .default = paste0("'", str_sub(x, 3,4)))
}
x_scale_global <- scale_x_continuous(limits=c(yearFirst-0.75, yearLast+0.75),
                                     breaks=seq(yearFirst, 2025, 5),
                                     labels = our_labels,
                                     minor_breaks=seq(1981, 2028, 1), 
                                     expand = c(0,0))
x_scale_marginal <- scale_x_continuous(limits=c(yearFirst-0.75, yearLast+0.75),
                                       breaks=seq(yearFirst, 2030, 5),
                                       labels = function(x) { paste0("'", str_sub(x, 3,4) ) },
                                       minor_breaks=NULL, 
                                       expand = c(0,0))

theme_marginal <- theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size=7.5),
        plot.caption = element_text(family="mono", size=7),
        plot.title.position = "plot",
        plot.subtitle = element_text(size=8.5) )

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    breaks <- unique(breaks) # Ensure unique breaks
    breaks[breaks == floor(breaks)] # Keep only integer values
  }
  return(fxn)
}


###########################################################
###########################################################
##  DO NOT EDIT THE CODE BELOW
##
##  It should all be automated given the variables
##    declared above
###########################################################



################################################################################
################################################################################
#####
################################################################################
################################################################################
## Plot of Global Findings
################################################################################
################################################################################

global_basin <- ibtracs_non_shorty_for_props %>%
  group_by(SEASON) %>%
  summarize(Hurricanes = sum(Hurricanes), 
            Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms) ) %>%
  ungroup() |>
  mutate(Prop_Major = Major_Storms/Total_Storms,
         Prop_Intense = Intense_Storms/Total_Storms)

global_basin_tall <- global_basin |>
  pivot_longer(-SEASON, names_to="Measure", values_to="Value") |>
  mutate(Measure = factor(Measure, 
                          levels=c("Total_Storms", "Hurricanes", "Major_Storms", "Prop_Major", "Intense_Storms", "Prop_Intense"),
                          labels=c("Tropical Cyclones",
                                   "Hurricanes",
                                   "Major Cyclones",
                                   "Proportion of Major Cyclones",
                                   "Intense Cyclones",
                                   "Proportion of Intense Cyclones") ) )

global_key_results_plot <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Basin=="Global") |>
  slice(rep(1:n(), each=4)) |>
  group_by(Measure=Series) |>
  mutate(Group = paste0("seg", rep(c("1", "2"), each=2) ),
         Season = case_when( (row_number()==1) ~ yearFirst-0.4,
                             (row_number()==2) ~ Time+0.4,
                             (row_number()==3) ~ Time+0.6,
                             (row_number()==4) ~ yearLast+0.4),
         Season = ifelse( (adj > sig_level & row_number()==2), 
                          Season + 1, 
                          Season),
         Value = case_when( (row_number()<=2 & adj < sig_level) ~ Start,
                            (row_number() >2 & adj < sig_level) ~ End,
                            .default = Mean) ) |>
  mutate(Measure = factor(Measure, 
                          levels=c("Cyclone Freq", "Hurricanes", "Major Freq", "Prop Major", "Intense Freq", "Prop Intense"),
                          labels=c("Tropical Cyclones",
                                   "Hurricanes",
                                   "Major Cyclones",
                                   "Proportion of Major Cyclones",
                                   "Intense Cyclones",
                                   "Proportion of Intense Cyclones") ) ) |>
  mutate(Value = ifelse((Measure %in% c("Proportion of Major Cyclones", "Proportion of Intense Cyclones") & adj < sig_level),
                        Value/100,
                        Value))
global_key_results_plot


plot_global_findings <- ggplot(global_basin_tall, aes(x=SEASON, y=Value) ) +
  geom_line(color=data_line) + 
  geom_point(color=data_dots) +
  facet_wrap(~Measure, nrow=6, scales="free_y") +
  geom_line(data=global_key_results_plot, aes(x=Season, y=Value, group=Group),
            linewidth=1.2, color=seg_line) +
  labs(title="Global Tropical Cyclone Record (1980-2025)",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_global

plot_global_findings

ggsave(filename="./plots/result_ibtracs_global_findings.png",
       plot=plot_global_findings,
       width=6.5, height=7, bg="white")










################################################################################
################################################################################
#####
################################################################################
################################################################################
## Total Cyclones per Basin
################################################################################
################################################################################

storm_counts <- ibtracs_non_shorty_for_props



storm_basin_regimes <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Series=="Cyclone Freq",
                Basin !="Global") |>
  slice(rep(1:n(), each=2)) |>
  group_by(BASIN=Basin) |>
  mutate(Regime = paste0("seg", c("1", "2") ),
         Min_year = case_when( (row_number()==1) ~ yearFirst-0.4,
                               (row_number()==2) ~ Time+0.6),
         Max_year = case_when( (row_number()==1) ~ Time+0.4,
                               (row_number()==2) ~ yearLast+0.04),
         Value = case_when( (row_number()==1 & adj < sig_level) ~ Start,
                            (row_number()==2 & adj < sig_level) ~ End,
                            .default = Mean) ) |>
  mutate(Max_year = ifelse( (adj > sig_level & row_number()==1), 
                            Max_year + 1, 
                            Max_year) ) |>
  mutate(BASIN = factor(BASIN,
                        levels=c("North Atlantic", "Eastern North Pacific", "Western North Pacific",
                                 "Northern Indian", "Southern Indian", "Southern Pacific") ) )


plot_basin_totals <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Total_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=Value, yend=Value, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Cyclone Frequency by Basin",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_marginal +
  scale_y_continuous(breaks = integer_breaks())

ggsave(filename="./plots/result_ibtracs_basin_nonshortiesfindings.png", 
       plot=plot_basin_totals,
       width=6.5, height=4.5, bg="white")





################################################################################
################################################################################
#####
################################################################################
################################################################################
## Hurricanes per Basin
################################################################################
################################################################################


storm_basin_regimes <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Series=="Hurricanes",
                Basin !="Global") |>
  slice(rep(1:n(), each=2)) |>
  group_by(BASIN=Basin) |>
  mutate(Regime = paste0("seg", c("1", "2") ),
         Min_year = case_when( (row_number()==1) ~ yearFirst-0.4,
                               (row_number()==2) ~ Time+0.6),
         Max_year = case_when( (row_number()==1) ~ Time+0.4,
                               (row_number()==2) ~ yearLast+0.04),
         Value = case_when( (row_number()==1 & adj < sig_level) ~ Start,
                            (row_number()==2 & adj < sig_level) ~ End,
                            .default = Mean) ) |>
  mutate(Max_year = ifelse( (adj > sig_level & row_number()==1), 
                            Max_year + 1, 
                            Max_year) ) |>
  mutate(BASIN = factor(BASIN,
                        levels=c("North Atlantic", "Eastern North Pacific", "Western North Pacific",
                                 "Northern Indian", "Southern Indian", "Southern Pacific") ) )



plot_basin_hurricanes <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Hurricanes), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=Value, yend=Value, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Hurricanes (Category 1+) by Basin",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  scale_y_continuous(breaks = integer_breaks()) +
  x_scale_marginal

ggsave(filename="./plots/result_ibtracs_basin_hurricanes_findings.png", 
       plot=plot_basin_hurricanes,
       width=6.5, height=4.5, bg="white")




################################################################################
################################################################################
#####
################################################################################
################################################################################
## Major Cyclones per Basin
################################################################################
################################################################################


storm_basin_regimes <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Series=="Major Freq",
                Basin !="Global") |>
  slice(rep(1:n(), each=2)) |>
  group_by(BASIN=Basin) |>
  mutate(Regime = paste0("seg", c("1", "2") ),
         Min_year = case_when( (row_number()==1) ~ yearFirst-0.4,
                               (row_number()==2) ~ Time+0.6),
         Max_year = case_when( (row_number()==1) ~ Time+0.4,
                               (row_number()==2) ~ yearLast+0.04),
         Value = case_when( (row_number()==1 & adj < sig_level) ~ Start,
                            (row_number()==2 & adj < sig_level) ~ End,
                            .default = Mean) ) |>
  mutate(Max_year = ifelse( (adj > sig_level & row_number()==1), 
                            Max_year + 1, 
                            Max_year) ) |>
  mutate(BASIN = factor(BASIN,
                        levels=c("North Atlantic", "Eastern North Pacific", "Western North Pacific",
                                 "Northern Indian", "Southern Indian", "Southern Pacific") ) )



plot_basin_major <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Major_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=Value, yend=Value, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Major Cyclones (Category 3+) by Basin",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  scale_y_continuous(breaks = integer_breaks()) +
  x_scale_marginal

ggsave(filename="./plots/result_ibtracs_basin_major_findings.png", 
       plot=plot_basin_major,
       width=6.5, height=4.5, bg="white")



################################################################################
################################################################################
#####
################################################################################
################################################################################
## Proportion Major Cyclones per Basin
################################################################################
################################################################################

storm_counts <- ibtracs_non_shorty_for_props |>
  mutate(Prop_Storms = Major_Storms/Total_Storms)


storm_basin_regimes <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Series=="Prop Major",
                Basin !="Global") |>
  slice(rep(1:n(), each=2)) |>
  group_by(BASIN=Basin) |>
  mutate(Regime = paste0("seg", c("1", "2") ),
         Min_year = case_when( (row_number()==1) ~ yearFirst-0.4,
                               (row_number()==2) ~ Time+0.6),
         Max_year = case_when( (row_number()==1) ~ Time+0.4,
                               (row_number()==2) ~ yearLast+0.04),
         Value = case_when( (row_number()==1 & adj < sig_level) ~ Start,
                            (row_number()==2 & adj < sig_level) ~ End,
                            .default = Mean),
         Value = ifelse(adj < sig_level, Value/100, Value) ) |>
  mutate(Max_year = ifelse( (adj > sig_level & row_number()==1), 
                            Max_year + 1, 
                            Max_year) ) |>
  mutate(BASIN = factor(BASIN,
                        levels=c("North Atlantic", "Eastern North Pacific", "Western North Pacific",
                                 "Northern Indian", "Southern Indian", "Southern Pacific") ) )




plot_basin_prop_major <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Prop_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=Value, yend=Value, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Proportion of Major Cyclones by Basin",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_marginal

ggsave(filename="./plots/result_ibtracs_basin_prop_major_findings.png", 
       plot=plot_basin_prop_major,
       width=6.5, height=4.5, bg="white")






################################################################################
################################################################################
#####
################################################################################
################################################################################
## Intense Cyclones per Basin
################################################################################
################################################################################


storm_basin_regimes <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Series=="Intense Freq",
                Basin !="Global") |>
  slice(rep(1:n(), each=2)) |>
  group_by(BASIN=Basin) |>
  mutate(Regime = paste0("seg", c("1", "2") ),
         Min_year = case_when( (row_number()==1) ~ yearFirst-0.4,
                               (row_number()==2) ~ Time+0.6),
         Max_year = case_when( (row_number()==1) ~ Time+0.4,
                               (row_number()==2) ~ yearLast+0.04),
         Value = case_when( (row_number()==1 & adj < sig_level) ~ Start,
                            (row_number()==2 & adj < sig_level) ~ End,
                            .default = Mean) ) |>
  mutate(Max_year = ifelse( (adj > sig_level & row_number()==1), 
                            Max_year + 1, 
                            Max_year) ) |>
  mutate(BASIN = factor(BASIN,
                        levels=c("North Atlantic", "Eastern North Pacific", "Western North Pacific",
                                 "Northern Indian", "Southern Indian", "Southern Pacific") ) )



plot_basin_intense <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Intense_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=Value, yend=Value, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Intense Cyclones (Category 4+) by Basin",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  scale_y_continuous(breaks = integer_breaks()) +
  theme_marginal +
  x_scale_marginal

ggsave(filename="./plots/result_ibtracs_basin_intense_findings.png", 
       plot=plot_basin_intense,
       width=6.5, height=4.5, bg="white")



################################################################################
################################################################################
#####
################################################################################
################################################################################
## Proportion Intense Cyclones per Basin
################################################################################
################################################################################

storm_counts <- ibtracs_non_shorty_for_props |>
  mutate(Prop_Storms = Intense_Storms/Total_Storms)

storm_basin_regimes <- ibtracs_key_results |>
  dplyr::filter(Method=="Changepoint") |>
  mutate(adj = p.adjust(p.value, method="fdr")) |>
  dplyr::filter(Series=="Prop Intense",
                Basin !="Global") |>
  slice(rep(1:n(), each=2)) |>
  group_by(BASIN=Basin) |>
  mutate(Regime = paste0("seg", c("1", "2") ),
         Min_year = case_when( (row_number()==1) ~ yearFirst-0.4,
                               (row_number()==2) ~ Time+0.6),
         Max_year = case_when( (row_number()==1) ~ Time+0.4,
                               (row_number()==2) ~ yearLast+0.04),
         Value = case_when( (row_number()==1 & adj < sig_level) ~ Start,
                            (row_number()==2 & adj < sig_level) ~ End,
                            .default = Mean),
         Value = ifelse(adj < sig_level, Value/100, Value) ) |>
  mutate(Max_year = ifelse( (adj > sig_level & row_number()==1), 
                            Max_year + 1, 
                            Max_year) ) |>
  mutate(BASIN = factor(BASIN,
                        levels=c("North Atlantic", "Eastern North Pacific", "Western North Pacific",
                                 "Northern Indian", "Southern Indian", "Southern Pacific") ) )



plot_basin_prop_intense <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Prop_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=Value, yend=Value, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Proportion of Intense Cyclones by Basin",
       subtitle=subtitle_text,
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_marginal

ggsave(filename="./plots/result_ibtracs_basin_prop_intense_findings.png", 
       plot=plot_basin_prop_intense,
       width=6.5, height=4.5, bg="white")



