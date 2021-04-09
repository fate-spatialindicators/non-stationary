## Fit Stationary and Non-Stationary Models to zero and positive responses ##

# new data,
# smaller knots, 600,
# delta lognormal instead of gamma
#remotes::install_github("pbs-assess/sdmTMB", ref = "epsilon-heterogeneity")
library(sdmTMB)
library(dplyr)

n_cutoff = 15 # generates ~ 600 knots
# Species of interest
species = read.csv("survey_data/species_list.csv")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

dat = readRDS("survey_data/all_data.rds")

# filter out a single species -- all we need is temp data
sub = dplyr::filter(dat,
                    scientific_name == unique(dat$scientific_name)[1])

sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2

# create calendar day
sub$month = as.numeric(substr(sub$date_yyyymmdd, 5, 6))
sub$day = as.numeric(substr(sub$date_yyyymmdd, 7, 8))
sub$Date = lubridate::parse_date_time(paste(sub$year, sub$month, sub$day), orders = "ymd")
sub$yday = lubridate::yday(sub$Date)

sub$calday = as.numeric(scale(sub$yday))
sub$calday2 = sub$calday^2

sub = dplyr::filter(sub, !is.na(temperature_at_gear_c_der))
spde = make_mesh(sub, c("lon", "lat"), cutoff = n_cutoff)

# total cpue
fit = try(sdmTMB(temperature_at_gear_c_der ~ 0 + depth_scaled + depth_scaled2 +
                   calday + calday2 + as.factor(year),
                    data = sub, time = "year", spde = spde)
)

# scale time
sub$time = sub$year - min(sub$year) + 1
fit_ll = try(sdmTMB(temperature_at_gear_c_der ~ 0 + depth_scaled + depth_scaled2 +
                      calday + calday2 + as.factor(year),
                       data = sub, time = "year", spde = spde, epsilon_predictor = "time")
)

temp = dplyr::group_by(sub, year) %>%
  dplyr::summarize(mean_temp = mean(temperature_at_gear_c_der,na.rm=T))
temp$mean_temp = scale(temp$mean_temp)
sub = dplyr::left_join(sub, temp)
fit_ll_temp = try(sdmTMB(temperature_at_gear_c_der ~ 0 + depth_scaled + depth_scaled2 +
                           calday + calday2 + as.factor(year),
                    data = sub, time = "mean_temp", spde = spde, epsilon_predictor = "time")
)
# drop out high memory objects
fit$data = NULL
fit$tmb_data = NULL
fit$spde = NULL
fit_ll$data = NULL
fit_ll$tmb_data = NULL
fit_ll$spde = NULL

if (!dir.exists("output")) {dir.create("output")}
save.image(file=paste0("output/temp_all_models.RData"))
