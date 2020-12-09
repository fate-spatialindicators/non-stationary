## Fit Non-Stationary Models ##
# [sdmTMB: non-stationary/log-linear model (process variance)]

#install.packages("remotes")
#remotes::install_github("pbs-assess/sdmTMB", ref = "epsilon-heterogeneity")

library(sdmTMB)
library(sp)
library(ggplot2)
library(ggsidekick)
library(dplyr)

# data
dat = readRDS("data/Darkblotched rockfish_expanded.rds")

# [index common name]

# make plots of variability in residuals for models of positive catch rates

# juveniles
sub_juv = dplyr::filter(dat, juv_cpue_kg_km2 > 0)
sub_juv$depth_scaled = as.numeric(scale(sub_juv$depth_m))
sub_juv$depth_scaled2 = sub_juv$depth_scaled^2
spde = make_mesh(sub_juv, c("lon", "lat"), cutoff = 10)
# check knots with spde$mesh$n
pos_juv_fit = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = sub_juv, time = "year", spde = spde, family = Gamma(link = "log"))
sub_juv$resid = sub_juv$juv_cpue_kg_km2 - predict(pos_juv_fit)$est
# plot of juvenile residual variation over time
g1 = dplyr::group_by(sub_juv, year) %>%
  dplyr::summarize(sd = sd(resid)) %>%
  ggplot(aes(year, log(sd))) + geom_point(size=2) + geom_line() + xlab("Year") + ylab("ln(Std. Dev.) Residuals Juveniles CPUE") +
  scale_x_continuous(breaks = seq(from = 2003, to = 2018, by = 1)) + ylim(-2.5,2.5) + theme_sleek()

#pdf("plots/Darkblotched rockfish_juvs_residuals_stationary.pdf")
g1
#dev.off()

# adults
sub_ad = dplyr::filter(dat, adult_cpue_kg_km2 > 0)
sub_ad$depth_scaled = as.numeric(scale(sub_ad$depth_m))
sub_ad$depth_scaled2 = sub_ad$depth_scaled^2
spde = make_mesh(sub_ad, c("lon", "lat"), cutoff = 10)
pos_ad_fit = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                    data = sub_ad, time = "year", spde = spde, family = Gamma(link = "log"))
sub_ad$resid = sub_ad$adult_cpue_kg_km2 - predict(pos_ad_fit)$est
# plot of adult residual variation over time
g2 = dplyr::group_by(sub_ad, year) %>%
  dplyr::summarize(sd = sd(resid)) %>%
  ggplot(aes(year, log(sd))) + geom_point(size=2) +geom_line() + xlab("Year") + ylab("ln(Std. Dev.) Residuals Adults CPUE") +
  scale_x_continuous(breaks = seq(from = 2003, to = 2018, by = 1)) + ylim(6.5,10.5) + theme_sleek()

#pdf("plots/Darkblotched rockfish_adults_residuals_stationary.pdf")
g2
#dev.off()

# fit non-stationary Tweedie model to zero and non-zero catches CPUE
dat$depth_scaled = as.numeric(scale(dat$depth_m))
dat$depth_scaled2 = dat$depth_scaled^2
spde = make_mesh(dat, c("lon", "lat"), cutoff = 10)

# adults and juveniles combined
total_fit = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                     epsilon_model = "loglinear")

dat$resid = dat$cpue_kg_km2 - predict(total_fit)$est 
# probably better to assign the what's returned by predict() to an object first unless this is the only operation you plan

# juveniles
juv_fit = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                   data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                   epsilon_model = "loglinear")

dat$resid = dat$juv_cpue_kg_km2 - predict(juv_fit)$est

# adults
ad_fit = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                 data = dat, time = "year", spde = spde, family = tweedie(link = "log"),
                 epsilon_model = "loglinear")

dat$resid = dat$adult_cpue_kg_km2 - predict(ad_fit)$est

# then you would do the same as the above 3 models for the stationary case....
# adults and juveniles combined
total_fit_stationary = sdmTMB(cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                   data = dat, time = "year", spde = spde, family = tweedie(link = "log"))

dat$resid = dat$cpue_kg_km2 - predict(total_fit_stationary)$est 

# juveniles
juv_fit_stationary = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                 data = dat, time = "year", spde = spde, family = tweedie(link = "log"))

dat$resid = dat$juv_cpue_kg_km2 - predict(juv_fit_stationary)$est

# adults
ad_fit_stationary = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                data = dat, time = "year", spde = spde, family = tweedie(link = "log"))

dat$resid = dat$adult_cpue_kg_km2 - predict(ad_fit_stationary)$est
