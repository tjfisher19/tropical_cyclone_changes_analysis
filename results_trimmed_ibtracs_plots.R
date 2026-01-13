################################################
##
##  results_trimmed_plots.R
##
##  This file builds the plots of the change
##    point results. If the found p-value
##    is less than 1/3 (the "likely" finding)
##    we plot the suggested segmentation base
##    on the HURSAT dataset.
##
##  The first plot is all the findings at the
##    global level. Then the subsequent plots
##    are the results within each of the 6 basins
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-01-13


library(tidyverse)
library(patchwork)

rm(list=ls())
yearFirst <- 1981
yearLast <- 2017


#######################################
## Load data & Change point findings

load("./data/ibtracs_nonShortiesForProportions.RData")
load("./data/results_trimmed.RData")

#########################################
## PLOT Parameters and
##   customization of axes and such
##   

data_dots <- "gray50"
data_line <- "gray60"
seg_line <- "gray20"
caption_text <- "Source: International Best Track Archive for Climate Stewardship\nhttps://www.ncei.noaa.gov/products/international-best-track-archive"

non_shorty_for_props <- non_shorty_for_props |>
  dplyr::filter(SEASON >= yearFirst,
                SEASON <= yearLast)

our_labels <- function(x) {
  case_when(x==yearFirst ~ paste0(yearFirst),
            x==2017 ~ "2017",
            x==2000 ~ "2000",
            .default = paste0("'", str_sub(x, 3,4)))
}
x_scale_global <- scale_x_continuous(limits=c(yearFirst-0.75, yearLast+0.75),
                                     breaks=seq(yearFirst, 2017, 3),
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
        axis.text = element_text(size=8),
        plot.caption = element_text(family="mono"),
        plot.title.position = "plot")

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    breaks <- unique(breaks) # Ensure unique breaks
    breaks[breaks == floor(breaks)] # Keep only integer values
  }
  return(fxn)
}

################################################################################
################################################################################
#####
################################################################################
################################################################################
## Plot of Global Findings
################################################################################
################################################################################

global_basin <- non_shorty_for_props %>%
  group_by(SEASON) %>%
  summarize(Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms) ) %>%
  ungroup() |>
  mutate(Prop_Major = Major_Storms/Total_Storms,
         Prop_Intense = Intense_Storms/Total_Storms)

global_basin_tall <- global_basin |>
  pivot_longer(-SEASON, names_to="Measure", values_to="Value") |>
  mutate(Measure = factor(Measure, 
                          levels=c("Total_Storms", "Major_Storms", "Prop_Major", "Intense_Storms", "Prop_Intense"),
                          labels=c("Total Non-Shorties",
                                   "Major Storms",
                                   "Proportion of Major Storms",
                                   "Intense Storms",
                                   "Proportion of Intense Storms") ) )

# Total Storms - change at time 20
trimmed_global_results$pois_chpt_test
# Major Storms - no change
trimmed_global_results$major_pois_chpt_test
# Proportion Major Storms - change at time 21
trimmed_global_results$prop_major_chpt_test
# Intense Storms - no change
trimmed_global_results$intense_pois_chpt_test
# Proportion Intense Storms - change at 21
trimmed_global_results$prop_intense_chpt_test

mean1 <- mean(global_basin$Total_Storms[1:trimmed_global_results$pois_chpt_test[2]])
mean2 <- mean(global_basin$Total_Storms[(trimmed_global_results$pois_chpt_test[2]+1):dim(global_basin)[1]])
df_total_seg <- data.frame(
  SEASON = c(min(global_basin$SEASON)-0.4, global_basin$SEASON[trimmed_global_results$pois_chpt_test[2]]+0.4,
             global_basin$SEASON[trimmed_global_results$pois_chpt_test[2]]+0.6, max(global_basin$SEASON)+0.4),
  Group = c(rep("seg1", 2),
            rep("seg2", 2) ),
  Mean = c(rep(mean1, 2),
           rep(mean2, 2) ),
  Measure = "Total_Storms" )

df_major_seg <- data.frame(
  SEASON = c(min(global_basin$SEASON)-0.4, max(global_basin$SEASON)+0.4),
  Group = c(rep("seg2", 2) ),
  Mean = c(rep(mean(global_basin$Major_Storms), 2) ),
  Measure = "Major_Storms")


mean1 <- mean(global_basin$Intense_Storms[1:trimmed_global_results$intense_pois_chpt_test[2]])
mean2 <- mean(global_basin$Intense_Storms[(trimmed_global_results$intense_pois_chpt_test[2]+1):dim(global_basin)[1]])
df_intense_seg <- data.frame(
  SEASON = c(min(global_basin$SEASON)-0.4, global_basin$SEASON[trimmed_global_results$intense_pois_chpt_test[2]]+0.4,
             global_basin$SEASON[trimmed_global_results$intense_pois_chpt_test[2]]+0.6, max(global_basin$SEASON)+0.4),
  Group = c(rep("seg1", 2),
            rep("seg2", 2) ),
  Mean = c(rep(mean1, 2),
           rep(mean2, 2) ),
  Measure = "Intense_Storms"
)


mean1 <- sum(global_basin$Major_Storms[1:trimmed_global_results$prop_major_chpt_test[2]])/
  sum(global_basin$Total_Storms[1:trimmed_global_results$prop_major_chpt_test[2]])
mean2 <- sum(global_basin$Major_Storms[(trimmed_global_results$prop_major_chpt_test[2]+1):dim(global_basin)[1]])/
  sum(global_basin$Total_Storms[(trimmed_global_results$prop_major_chpt_test[2]+1):dim(global_basin)[1]])
df_prop_major_seg <- data.frame(
  SEASON = c(min(global_basin$SEASON)-0.4, global_basin$SEASON[trimmed_global_results$prop_major_chpt_test[2]]+0.4,
             global_basin$SEASON[trimmed_global_results$prop_major_chpt_test[2]]+0.6, max(global_basin$SEASON)+0.4),
  Group = c(rep("seg1", 2),
            rep("seg2", 2) ),
  Mean = c(rep(mean1, 2),
           rep(mean2, 2) ),
  Measure = "Prop_Major"
)



mean1 <- sum(global_basin$Intense_Storms[1:trimmed_global_results$prop_intense_chpt_test[2]])/
  sum(global_basin$Total_Storms[1:trimmed_global_results$prop_intense_chpt_test[2]])
mean2 <- sum(global_basin$Intense_Storms[(trimmed_global_results$prop_intense_chpt_test[2]+1):dim(global_basin)[1]])/
  sum(global_basin$Total_Storms[(trimmed_global_results$prop_intense_chpt_test[2]+1):dim(global_basin)[1]])
df_prop_intense_seg <- data.frame(
  SEASON = c(min(global_basin$SEASON)-0.4, global_basin$SEASON[trimmed_global_results$prop_intense_chpt_test[2]]+0.4,
             global_basin$SEASON[trimmed_global_results$prop_intense_chpt_test[2]]+0.6, max(global_basin$SEASON)+0.4),
  Group = c(rep("seg1", 2),
            rep("seg2", 2) ),
  Mean = c(rep(mean1, 2),
           rep(mean2, 2) ),
  Measure = "Prop_Intense"
)

df_segments <- bind_rows(
  df_total_seg, df_major_seg, df_prop_major_seg, df_intense_seg, df_prop_intense_seg) |>
  mutate(Measure = factor(Measure, levels=c("Total_Storms", "Major_Storms", "Prop_Major", "Intense_Storms", "Prop_Intense"),
                          labels=c("Total Non-Shorties",
                                   "Major Storms",
                                   "Proportion of Major Storms",
                                   "Intense Storms",
                                   "Proportion of Intense Storms")))

plot_global_findings <- ggplot(global_basin_tall, aes(x=SEASON, y=Value) ) +
  geom_line(color=data_line) + 
  geom_point(color=data_dots) +
  facet_wrap(~Measure, nrow=5, scales="free_y") +
  geom_line(data=df_segments, aes(x=SEASON, y=Mean, group=Group),
            linewidth=1.2, color=seg_line) +
  labs(title="Global Tropical Cyclone Record (1981-2017)",
       subtitle="Segmentations based on SCUSUM Change Point Test",
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_global

plot_global_findings

ggsave(filename="./plots/result_trimmed_global_findings.png",
       plot=plot_global_findings,
       width=6.5, height=7, bg="white")










###########################################################
###########################################################
##  DO NOT EDIT THE CODE BELOW
##
##  It should all be automated given the variables
##    declared above
###########################################################

####################################################
######### Proportions
####################################################

title_text <- "Proportion of Major Storms to all non-'shorty' Cyclones per Season"
subtitle_text <- "Proportion of storms attaining at least 96 knots, with BIC segmentation restricted to regimes of at least"
subtitle_text_norestrict <- "Proportion of storms attaining at least 96 knots with BIC segmentation"
file_name <- "propMajor"

N <- length(unique(non_shorty_for_props$SEASON))
chpts_prop <- 1:(N-1)
chpts_lambda <- N:(2*(N-1))

global_counts <- non_shorty_for_props %>%
  group_by(SEASON) %>%
  summarize(Major_Storms = sum(Major_Storms),
            Intense_Storms = sum(Intense_Storms),
            Total_Storms = sum(Total_Storms) ) %>%
  ungroup() |>
  mutate(Prop_Major = Major_Storms/Total_Storms,
         Prop_Intense = Intense_Storms/Total_Storms)


################################################################################
################################################################################
#####
################################################################################
################################################################################
## Total Storms per Basin
################################################################################
################################################################################

trimmed_na_results$pois_chpt_test
trimmed_ep_results$pois_chpt_test
trimmed_wp_results$pois_chpt_test
trimmed_ni_results$pois_chpt_test
trimmed_si_results$pois_chpt_test
trimmed_sp_results$pois_chpt_test

storm_counts <- non_shorty_for_props |>
  mutate(Prop_Storms = Major_Storms/Total_Storms)

storm_basin_regimes <- storm_counts %>%
  ungroup() %>%
  arrange(BASIN, SEASON) %>%
  mutate(Regime = c(findInterval(1:N, trimmed_na_results$pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_ep_results$pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_wp_results$pois_chpt_test[2]+1),
                    findInterval(1:N, NULL),
                    findInterval(1:N, trimmed_si_results$pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_sp_results$pois_chpt_test[2]+1)) ) %>%
  group_by(BASIN, Regime) %>%
  summarize(Regime_Length = max(SEASON)-min(SEASON)+1,
            AvgTotalStorms = mean(Total_Storms),
            Min_year = min(SEASON)-0.4,
            Max_year = max(SEASON)+0.4 )


plot_basin_totals <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Total_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=AvgTotalStorms, yend=AvgTotalStorms, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Non-Shorties by Basin",
       subtitle="Segmentations based on SCUSUM Change Point Test",
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_marginal +
  scale_y_continuous(breaks = integer_breaks())

ggsave(filename="./plots/result_trimmed_basin_nonshortiesfindings.png", 
       plot=plot_basin_totals,
       width=6.5, height=4.5, bg="white")





################################################################################
################################################################################
#####
################################################################################
################################################################################
## Major Storms per Basin
################################################################################
################################################################################

trimmed_na_results$major_pois_chpt_test
trimmed_ep_results$major_pois_chpt_test
trimmed_wp_results$major_pois_chpt_test
trimmed_ni_results$major_pois_chpt_test
trimmed_si_results$major_pois_chpt_test
trimmed_sp_results$major_pois_chpt_test

storm_counts <- non_shorty_for_props |>
  mutate(Prop_Storms = Major_Storms/Total_Storms)

storm_basin_regimes <- storm_counts %>%
  ungroup() %>%
  arrange(BASIN, SEASON) %>%
  mutate(Regime = c(findInterval(1:N, trimmed_na_results$major_pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_ep_results$major_pois_chpt_test[2]+1),
                    findInterval(1:N, NULL),
                    findInterval(1:N, NULL),
                    findInterval(1:N, NULL),
                    findInterval(1:N, NULL) ) ) %>%
  group_by(BASIN, Regime) %>%
  summarize(Regime_Length = max(SEASON)-min(SEASON)+1,
            AvgMajorStorms = mean(Major_Storms),
            Min_year = min(SEASON)-0.4,
            Max_year = max(SEASON)+0.4 )


plot_basin_major <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Major_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=AvgMajorStorms, yend=AvgMajorStorms, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Major Cyclones (Category 3+) by Basin",
       subtitle="Segmentations based on SCUSUM Change Point Test",
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  scale_y_continuous(breaks = integer_breaks()) +
  x_scale_marginal

ggsave(filename="./plots/result_trimmed_basin_major_findings.png", 
       plot=plot_basin_major,
       width=6.5, height=4.5, bg="white")



################################################################################
################################################################################
#####
################################################################################
################################################################################
## Proportion Major Storms per Basin
################################################################################
################################################################################

trimmed_na_results$prop_major_chpt_test
trimmed_ep_results$prop_major_chpt_test
trimmed_wp_results$prop_major_chpt_test
trimmed_ni_results$prop_major_chpt_test
trimmed_si_results$prop_major_chpt_test
trimmed_sp_results$prop_major_chpt_test

storm_counts <- non_shorty_for_props |>
  mutate(Prop_Storms = Major_Storms/Total_Storms)

storm_basin_regimes <- storm_counts %>%
  ungroup() %>%
  arrange(BASIN, SEASON) %>%
  mutate(Regime = c(findInterval(1:N, NULL),
                    findInterval(1:N, NULL),
                    findInterval(1:N, trimmed_wp_results$prop_major_chpt_test[2]+1),
                    findInterval(1:N, NULL),
                    findInterval(1:N, trimmed_si_results$prop_major_chpt_test[2]+1),
                    findInterval(1:N, trimmed_sp_results$prop_major_chpt_test[2]+1) ) ) |>
  group_by(BASIN, Regime) %>%
  summarize(Regime_Length = max(SEASON)-min(SEASON)+1,
            PropMajorStorms = sum(Major_Storms)/sum(Total_Storms),
            Min_year = min(SEASON)-0.4,
            Max_year = max(SEASON)+0.4 )


plot_basin_prop_major <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Prop_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=PropMajorStorms, yend=PropMajorStorms, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Proportion of Major Cyclones by Basin",
       subtitle="Segmentations based on SCUSUM Change Point Test",
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_marginal

ggsave(filename="./plots/result_trimmed_basin_prop_major_findings.png", 
       plot=plot_basin_prop_major,
       width=6.5, height=4.5, bg="white")






################################################################################
################################################################################
#####
################################################################################
################################################################################
## Intense Storms per Basin
################################################################################
################################################################################

trimmed_na_results$intense_pois_chpt_test
trimmed_ep_results$intense_pois_chpt_test
trimmed_wp_results$intense_pois_chpt_test
trimmed_ni_results$intense_pois_chpt_test
trimmed_si_results$intense_pois_chpt_test
trimmed_sp_results$intense_pois_chpt_test

storm_counts <- non_shorty_for_props |>
  mutate(Prop_Storms = Intense_Storms/Total_Storms)

storm_basin_regimes <- storm_counts %>%
  ungroup() %>%
  arrange(BASIN, SEASON) %>%
  mutate(Regime = c(findInterval(1:N, trimmed_na_results$intense_pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_ep_results$intense_pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_wp_results$intense_pois_chpt_test[2]+1),
                    findInterval(1:N, NULL),
                    findInterval(1:N, trimmed_si_results$intense_pois_chpt_test[2]+1),
                    findInterval(1:N, trimmed_sp_results$intense_pois_chpt_test[2]+1) ) ) %>%
  group_by(BASIN, Regime) %>%
  summarize(Regime_Length = max(SEASON)-min(SEASON)+1,
            AvgIntenseStorms = mean(Intense_Storms),
            Min_year = min(SEASON)-0.4,
            Max_year = max(SEASON)+0.4 )


plot_basin_intense <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Intense_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=AvgIntenseStorms, yend=AvgIntenseStorms, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Intense Cyclones (Category 4+) by Basin",
       subtitle="Segmentations based on SCUSUM Change Point Test",
       caption=caption_text) +
  theme_minimal() + 
  scale_y_continuous(breaks = integer_breaks()) +
  theme_marginal +
  x_scale_marginal

ggsave(filename="./plots/result_trimmed_basin_intense_findings.png", 
       plot=plot_basin_intense,
       width=6.5, height=4.5, bg="white")



################################################################################
################################################################################
#####
################################################################################
################################################################################
## Proportion Intense Storms per Basin
################################################################################
################################################################################

trimmed_na_results$prop_intense_chpt_test
trimmed_ep_results$prop_intense_chpt_test
trimmed_wp_results$prop_intense_chpt_test
trimmed_ni_results$prop_intense_chpt_test
trimmed_si_results$prop_intense_chpt_test
trimmed_sp_results$prop_intense_chpt_test

storm_counts <- non_shorty_for_props |>
  mutate(Prop_Storms = Intense_Storms/Total_Storms)

storm_basin_regimes <- storm_counts %>%
  ungroup() %>%
  arrange(BASIN, SEASON) %>%
  mutate(Regime = c(findInterval(1:N, NULL),
                    findInterval(1:N, trimmed_ep_results$prop_intense_chpt_test[2]+1),
                    findInterval(1:N, trimmed_wp_results$prop_intense_chpt_test[2]+1),
                    findInterval(1:N, NULL),
                    findInterval(1:N, trimmed_si_results$prop_intense_chpt_test[2]+1),
                    findInterval(1:N, trimmed_sp_results$prop_intense_chpt_test[2]+1) ) ) |>
  group_by(BASIN, Regime) %>%
  summarize(Regime_Length = max(SEASON)-min(SEASON)+1,
            PropIntenseStorms = sum(Intense_Storms)/sum(Total_Storms),
            Min_year = min(SEASON)-0.4,
            Max_year = max(SEASON)+0.4 )


plot_basin_prop_intense <- ggplot(storm_counts) + 
  geom_line(aes(x=SEASON, y=Prop_Storms), color=data_line) + 
  geom_segment(data=storm_basin_regimes, aes(x=Min_year, xend=Max_year,
                                             y=PropIntenseStorms, yend=PropIntenseStorms, group=Regime),
               color=seg_line, linewidth=1.25) +
  facet_wrap(~BASIN, scales="free") +
  labs(title="Changes in Proportion of Intense Cyclones by Basin",
       subtitle="Segmentations based on SCUSUM Change Point Test",
       caption=caption_text) +
  theme_minimal() + 
  theme_marginal +
  x_scale_marginal

ggsave(filename="./plots/result_trimmed_basin_prop_intense_findings.png", 
       plot=plot_basin_prop_intense,
       width=6.5, height=4.5, bg="white")



