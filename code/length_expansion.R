library(sdmTMB)
library(sp)
library(ggplot2)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)

species = read.csv("survey_data/species_list.csv") #data.frame(Scientific.Name = c("Anoplopoma fimbria","Sebastes alutus"))
bio = readRDS("survey_data/wcbts_bio_2019-08-01.rds")
haul = readRDS("survey_data/wcbts_haul_2019-08-01.rds")
catch = readRDS("survey_data/wcbts_catch_2019-08-01.rds")
names(catch) = tolower(names(catch))
names(bio) = tolower(names(bio))
names(haul) = tolower(names(haul))

bio$trawl_id = as.character(bio$trawl_id)
haul$trawl_id = as.character(haul$trawl_id)
haul$date_yyyymmdd = as.numeric(haul$date_yyyymmdd)
haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)



# join in bio and haul data
dat = dplyr::left_join(catch[,c("trawl_id","scientific_name","year","subsample_count",
  "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")], haul) %>%
  dplyr::left_join(bio) %>%
  filter(!is.na(length_cm), performance == "Satisfactory")

# filter species from list
dat = dplyr::filter(dat, scientific_name %in% species$Scientific.Name)
# bring in catch data which has total kg, filter haul performance and missing lengths
#dat = dat %>% left_join(catch[,c("trawl_id","scientific_name","year","subsample_count",
#  "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")]) %>%


# do spatial conversion
coordinates(dat) <- c("longitude_dd", "latitude_dd")
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
dat <- spTransform(dat, CRS(newproj))
dat = as.data.frame(dat)
dat$lon = dat$longitude_dd/1000
dat$lat = dat$latitude_dd/1000

# filter out species of interest
dat = dplyr::filter(dat, scientific_name == "Sebastes crameri")

# first-stage expansion from subsample to total sample biomass:
# for each trawl_id, (1) figure out weight of fish in juvenile length bin
# (2) if subsample weight == total weight, expansion factor = 1. (3) if
# subsample weight < total weight, expansion factor = total / subsample wt

# fit length-weight regression by year and sex to predict fish weights that have lengths only.
# note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)
fitted = dat %>%
  select(common_name, scientific_name, year, trawl_id, lon, lat,
         depth_m, o2_at_gear_ml_per_l_der, salinity_at_gear_psu_der, temperature_at_gear_c_der,
         subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
         individual_tracking_id, sex, length_cm, weight_kg) %>%
  group_nest(year, sex)  %>%
  mutate(
    model = map(data, ~ lm(log(weight_kg) ~ log(length_cm), data = .x)),
    tidied = map(model, tidy),
    augmented = map(model, augment),
    predictions = map2(data, model, modelr::add_predictions)
  )

# replace missing weights with predicted weights
dat = fitted %>%
  unnest(predictions) %>%
  select(-data, -model, -tidied, -augmented) %>%
  mutate(weight = ifelse(is.na(weight_kg), exp(pred), weight_kg))

# define length cutoff to define ontogenetic classes
juv_threshold = 15 # sablefish specific

# this just summarizes data at trawl_id level and sums up juv_weight
expanded = dplyr::group_by(dat, trawl_id) %>%
  dplyr::summarize(lon = lon[1], lat = lat[1], year = year[1],
    area_swept_ha_der = area_swept_ha_der[1],
    total_catch_wt_kg = total_catch_wt_kg[1],
    cpue_kg_km2 = cpue_kg_km2[1],
    subsample_wt_kg = subsample_wt_kg[1],
    juv_weight = sum(weight[which(length_cm < juv_threshold)]),
    adult_weight = sum(weight[which(length_cm > juv_threshold)])) %>%
  dplyr::filter(!is.na(total_catch_wt_kg), !is.na(area_swept_ha_der))
# expansion ratio is 1 for trawls with 100% subsampled catch. affects ~ 10% of trawls
expanded$ratio = 1
indx = which(expanded$subsample_wt_kg < expanded$total_catch_wt_kg)
expanded$ratio[indx] = expanded$total_catch_wt_kg[indx]/expanded$subsample_wt_kg[indx]



# calculate cpue for juvenile
expanded$juv_cpue = expanded$ratio * expanded$juv_weight / expanded$area_swept_ha_der


