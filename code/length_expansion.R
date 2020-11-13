library(sdmTMB)
library(sp)
library(ggplot2)
species = read.csv("survey_data/species_list.csv")
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
dat = dplyr::left_join(bio, haul)
# filter species from list
dat = dplyr::filter(dat, scientific_name %in% species$Scientific.Name)
# bring in catch data which has total kg
dat = dplyr::left_join(dat, catch[,c("trawl_id","scientific_name","year","subsample_count",
  "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")])

# filter out sablefish
dat = dplyr::filter(dat, scientific_name == "Anoplopoma fimbria")

# do spatial conversion
coordinates(dat) <- c("longitude_dd", "latitude_dd")
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
newproj = paste("+proj=utm +zone=10 ellps=WGS84")
dat <- spTransform(dat, CRS(newproj))
dat = as.data.frame(dat)
dat$lon = dat$longitude_dd/1000
dat$lat = dat$latitude_dd/1000

# expansion:
# for each trawlid, (1) figure out weight of fish in juvenile length bin
# (2) if subsample weight == total weight, expansion factor = 1. (3) if
# subssample weight < total weight, expansion factor = total / subsample wt

# filter out a few missing lengths
dat = dplyr::filter(dat, !is.na(length_cm))
# fit length-weight regression to predict fish weights that have lengths only
fit = lm(log(weight_kg) ~ log(length_cm), data=dat)
dat$pred_weight = exp(predict(fit, newdata=dat))
dat = dat %>%
  dplyr::mutate(weight = ifelse(is.na(length_cm), pred_weight, weight_kg))

juv_threshold = 29 # sablefish specific
# this just summarizes data at trawl_id level and sums up juv_weight
expanded = dplyr::group_by(dat, trawl_id) %>%
  dplyr::summarize(lon = lon[1], lat = lat[1], year = year[1],
    area_swept_ha_der = area_swept_ha_der[1],
    total_catch_wt_kg = total_catch_wt_kg[1],
    cpue_kg_km2 = cpue_kg_km2[1],
    subsample_wt_kg = subsample_wt_kg[1],
    juv_weight = sum(weight[which(length_cm < juv_threshold)])) %>%
  dplyr::filter(!is.na(total_catch_wt_kg))
# expansion ratio is 1 for trawls with 100% subsampled catch. affects ~ 10% of trawls
expanded$ratio = 1
indx = which(expanded$subsample_wt_kg < expanded$total_catch_wt_kg)
expanded$ratio[indx] = expanded$total_catch_wt_kg[indx]/expanded$subsample_wt_kg[indx]

# calculate cpue for juvenile
expanded$juv_cpue = expanded$ratio * expanded$juv_weight / expanded$area_swept_ha_der


