################################################################################
#### Project: SNF Field Experiment
#### Title:   Prepare soil data for Lena
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    30 March 2021
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
library(phyloseq)

## Code ----
scripts <- list.files(
  "./code",
  full.names = T
)

sapply(
  scripts,
  source
)

## Data ----
# MCC data
bac <- drake::readd(bacFull)$percGMPR
fun <- drake::readd(funFull)$percGMPR

# Home folder
tpFolder <- "/Users/tomwalker/Dropbox/projects/2018_transplant/SNF_field/"
soil2018 <- fread(paste0(tpFolder, "data/soil_data/2018_end/snf_field_2018_end_all_soil.csv"), data.table = F)
soil2019_1 <- fread(paste0(tpFolder, "data/soil_data/2019/snf_field_2019_early_all_soil.csv"), data.table = F)
soil2019_2 <- fread(paste0(tpFolder, "data/soil_data/2019/snf_field_2019_mid_all_soil.csv"), data.table = F)
soil2019_3 <- fread(paste0(tpFolder, "data/soil_data/2019/snf_field_2019_end_all_soil.csv"), data.table = F)


#### FORMAT MCC ----------------------------------------------------------------
## Select site ----
bacSub <- subset_samples(bac, site != "setup")
funSub <- subset_samples(fun, site != "setup")

## Do NMDS ----
bacNMDS <- do_nmds(bacSub)
funNMDS <- do_nmds(funSub)

## Extract scores etc ----
bacScores <- bacNMDS$points %>%
  as.data.frame %>%
  bind_cols(unphylo(bacSub@sam_data), .) %>%
  mutate(site_treat = paste(site, treat_a)) %>%
  select(-experiment)
funScores <- funNMDS$points %>%
  as.data.frame %>%
  bind_cols(unphylo(funSub@sam_data), .) %>%
  mutate(site_treat = paste(site, treat_a)) %>%
  select(-experiment)
allScores <- bind_cols(
  select(bacScores, site:rep_block),
  data.frame(
    fun_nmds1 = funNMDS$points[, 1],
    fun_nmds2 = funNMDS$points[, 2],
    bac_nmds1 = bacNMDS$points[, 1],
    bac_nmds2 = funNMDS$points[, 2]
  )
)
  
  
#### FORMAT SOIL ---------------------------------------------------------------

## Check order of high site (late 2019) ----
# original order
par(mfrow = c(2, 3))
plot(soil2019_1$cmic_ugC_g, soil2019_2$cmic_ugC_g)
plot(soil2019_3$cmic_ugC_g, soil2019_2$cmic_ugC_g)
plot(soil2019_3$cmic_ugC_g, soil2019_1$cmic_ugC_g)
# other order
plot(soil2019_1$cmic_ugC_g, soil2019_2$cmic_ugC_g)
plot(soil2019_3[c(11:60, 1:10), "cmic_ugC_g"], soil2019_2$cmic_ugC_g)
plot(soil2019_3[c(11:60, 1:10), "cmic_ugC_g"], soil2019_1$cmic_ugC_g)
# order is basically OK

## Average across three summer dates in 2019 ----
# average pool data
soil2019pools <- bind_rows(
  soil2019_1,
  select(soil2019_2, plot:no3_ugN_g),
  soil2019_3
) %>%
  group_by(plot, elevation, block, treatment, soil, plants) %>%
  summarise_at(vars(doc_ugC_g:no3_ugN_g), mean, na.rm = T)
# add enzyme data
soil2019 <- soil2019_2 %>%
  select(plot, GCS_nmol_h_g:APS_nmol_h_g) %>%
  left_join(soil2019pools, .)


#### EXPORT --------------------------------------------------------------------

write.table(
  x = allScores,
  file = "./exported/mcc_for_lena.csv",
  row.names = F, sep = ","
)
write.table(
  x = soil2019,
  file = "./exported/soilpools_2019_lena.csv",
  row.names = F, sep = ","
)
write.table(
  x = soil2018,
  file = "./exported/soilpools_2018_lena.csv",
  row.names = F, sep = ","
)

  
  
  
  
  
  
  


