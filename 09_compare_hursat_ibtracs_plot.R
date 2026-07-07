#######################################################
##
##  09_compare_hursat_ibtracs_plot.R
##
##  This code makes a gray-scale better formatted
##    density plot (see 03_ibtracs_vs_hursat report)
##    that is included in the supplemental material
##    of the paper
##
##  Author: Tom Fisher (fishert4@miamioh.edu)
##
##  Code tested on 2026-07-07



library(tidyverse)


load("data/ibtracs_storm_classifications.RData")
load("data/hursat_storm_classifications.RData")
ibtracs_storms_fct <- ibtracs_storms_fct |>
  filter(BASIN!="Southern Atlantic")
hursat_storms_fct <- hursat_storms_fct |>
  filter(BASIN!="Southern Atlantic")

compare_wind_df <- bind_rows(
  ibtracs_storms_fct |> 
    dplyr::filter(Duration_Class_Streak=="Non_Shorty") |>
    mutate(Source="IBTrACS"),
  
  hursat_storms_fct |>
    dplyr::filter(Duration_Class_Streak=="Non_Shorty") |>
    mutate(Source="ADT-HURSAT")
) |>
  mutate(Wind_ms = Wind/1.944)

plot_compare_wind_distro <- ggplot(compare_wind_df, aes(x=Wind_ms, fill=Source) ) +
  geom_density(alpha=0.4) +
  theme_minimal() +
  scale_fill_manual(name="Source:", values=c("black", "gray80") ) +
  labs(x="Maximum Windspeed (m/s)",
       title="Distribution of Maximum Windspeeds",
       subtitle="Includes All Storms at that reached Tropical Cyclone Strength for 3+ days",
       caption="Sources: ADT-HURSAT and IBTrACS\nhttps://www.ncei.noaa.gov/products/advanced-dvorak-technique-hurricane-satellite\nhttps://www.ncei.noaa.gov/products/international-best-track-archive"
  ) + 
  theme_minimal() +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size=7.5),
        plot.caption = element_text(family="mono", size=7),
        plot.title.position = "plot",
        plot.subtitle = element_text(size=8.5) )

ggsave(plot=plot_compare_wind_distro,
       filename="./plots/compare_ibtracs_hursat.png",
       width=6.5, height=4.25, bg="white")
