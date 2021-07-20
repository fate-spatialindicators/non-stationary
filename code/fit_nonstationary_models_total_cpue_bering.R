## Fit Stationary and Non-Stationary Models to Bering Cod data ##

#remotes::install_github("pbs-assess/sdmTMB@pc-prior")
library(sdmTMB)
library(dplyr)
library(sp)
library(rgdal)
if (!dir.exists("output")) dir.create("output")

n_cutoff <- 30

dat <- readRDS("survey_data/EBS_NBS_Index.RDS")
dat <- dat %>% dplyr::filter(SPECIES_CODE == 21740 & REGION == "EBS") # Pcod 21720; pollock 21740
dat <-  dplyr::transmute(dat,
                         cpue_kg_km2 = wCPUE*100,
                         year = YEAR,
                         lat = START_LATITUDE,
                         lon = START_LONGITUDE)

coordinates(dat) <- c("lon", "lat")
proj4string(dat) <- CRS("+proj=longlat +datum=WGS84")
dat <- spTransform(dat, CRS("+proj=utm +zone=3"))
dat <- as.data.frame(dat)
dat$lon <- dat$lon/1000 # scale units to km
dat$lat <- dat$lat/1000

dat$year <- as.factor(dat$year)
dat$time <- as.numeric(dat$year) - floor(mean(unique(as.numeric(dat$year))))

# AR1 spatiotemporal:
fit_models_ar1 <- function(sub) {
  spde <- make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff)
  ad_fit <- tryCatch({
    sdmTMB(cpue_kg_km2 ~ 0 + year,
      ar1_fields = TRUE, include_spatial = FALSE,
      data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
      matern_prior_E = c(50, 0.05, 1, 0.05),
      nlminb_loops = 2, newton_steps = 1
    )
  }, error = function(e) NA)
  ad_fit <- refit_model_if_needed(ad_fit)
  ad_fit_ll <- tryCatch({
    sdmTMB(cpue_kg_km2 ~ 0 + year,
      data = sub, time = "year", spde = spde,
      ar1_fields = TRUE, include_spatial = FALSE,
      family = tweedie(link = "log"), epsilon_predictor = "time",
      matern_prior_E = c(50, 0.05, 1, 0.05),
      nlminb_loops = 2, newton_steps = 1
    )}, error = function(e) NA)
  ad_fit_ll <- refit_model_if_needed(ad_fit_ll)
  save(ad_fit, ad_fit_ll,
    file = "output/pollock_bering_ar1.RData")
}

refit_model_if_needed <- function(m) {
  if (!is.na(m[[1]])) {
    if (max(m$gradients) > 0.01) {
      m <- tryCatch({
        sdmTMB::run_extra_optimization(m,
          nlminb_loops = 1L,
          newton_steps = 1L
        )
      }, error = function(e) m)
    }
  }
  m
}

fit_models_ar1(dat)
