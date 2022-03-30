################################################################################
#### Project: SNF Field Experiment
#### Title:   Analyse MCC in starting soils
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    26 March 2021
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())
# R session options (no factors, bias against scientific #s)
options(
stringsAsFactors = F,
scipen = 6
)

## Libraries ----
# standard library set
library(tidyverse)
library(data.table)

## Code ----
scripts <- list.files(
  "./code",
  full.names = T
)

sapply(
  scripts,
  source
)

## Load data ----
bac <- drake::readd(bacFull)$percGMPR
fun <- drake::readd(funFull)$percGMPR

## Select data ----
bacSub <- subset_samples(bac, site == "setup")
funSub <- subset_samples(fun, site == "setup")

#### NMDS ----------------------------------------------------------------------

## Do NMDS ----
bacNMDS <- do_nmds(bacSub)
funNMDS <- do_nmds(funSub)

## Extract scores etc ----
bacScores <- bacNMDS$points %>%
  as.data.frame %>%
  bind_cols(unphylo(bacSub@sam_data), .)
funScores <- funNMDS$points %>%
  as.data.frame %>%
  bind_cols(unphylo(funSub@sam_data), .)

## Get biplot data ----
bacPlotData <- make_biplot_data(
  PC1 = bacScores$MDS1,
  PC2 = bacScores$MDS2,
  groupVar = bacScores$treat_a
)
funPlotData <- make_biplot_data(
  PC1 = funScores$MDS1,
  PC2 = funScores$MDS2,
  groupVar = funScores$treat_a
)

## Plot ----
ggplot(bacPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  # scale_colour_manual(values = myCols) +
  # guides(col = "none") +
  aes(x = PC1, y = PC2, col = groupVar, xend = meanPC1, yend = meanPC2) +
  geom_segment() +
  geom_point() +
  geom_point(aes(x = meanPC1, y = meanPC2), shape = 21, fill = "white", size = 3)
ggplot(funPlotData) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  # scale_colour_manual(values = myCols) +
  # guides(col = "none") +
  aes(x = PC1, y = PC2, col = groupVar, xend = meanPC1, yend = meanPC2) +
  geom_segment() +
  geom_point() +
  geom_point(aes(x = meanPC1, y = meanPC2), shape = 21, fill = "white", size = 3)