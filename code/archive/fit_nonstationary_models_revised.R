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
sub = dplyr::filter(dat, juv_cpue_kg_km2 > 0)
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2
spde = make_mesh(sub, c("lon", "lat"), cutoff = 10)
# check knots with spde$mesh$n
pos_juv_fit = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = sub, time = "year", spde = spde, family = Gamma(link = "log"))
sub$resid = sub$juv_cpue_kg_km2 - predict(pos_juv_fit)$est
# plot of juvenile residual variation over time
g1 = dplyr::group_by(sub, year) %>%
  dplyr::summarize(sd = sd(resid)) %>%
  ggplot(aes(year, log(sd))) + geom_point(size=2) + geom_line() + xlab("Year") + ylab("ln(Std. Dev.) Residuals Juveniles CPUE") +
  scale_x_continuous(breaks = seq(from = 2003, to = 2018, by = 1)) + ylim(-2.5,2.5) + theme_sleek()

#pdf("plots/Darkblotched rockfish_juvs_residuals_stationary.pdf")
g1
#dev.off()

# adults
sub = dplyr::filter(dat, adult_cpue_kg_km2 > 0)
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2
spde = make_mesh(sub, c("lon", "lat"), cutoff = 10)
pos_ad_fit = sdmTMB(adult_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                    data = sub, time = "year", spde = spde, family = Gamma(link = "log"))
sub$resid = sub$adult_cpue_kg_km2 - predict(pos_ad_fit)$est
# plot of adult residual variation over time
g2 = dplyr::group_by(sub, year) %>%
  dplyr::summarize(sd = sd(resid)) %>%
  ggplot(aes(year, log(sd))) + geom_point(size=2) +geom_line() + xlab("Year") + ylab("ln(Std. Dev.) Residuals Adults CPUE") +
  scale_x_continuous(breaks = seq(from = 2003, to = 2018, by = 1)) + ylim(6.5,10.5) + theme_sleek()

#pdf("plots/Darkblotched rockfish_adults_residuals_stationary.pdf")
g2
#dev.off()

# fit non-stationary Tweedie model to total CPUE
sub$depth_scaled = as.numeric(scale(sub$depth_m))
sub$depth_scaled2 = sub$depth_scaled^2
spde = make_mesh(sub, c("lon", "lat"), cutoff = 10)
juv_tweedie = sdmTMB(juv_cpue_kg_km2 ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
                     data = sub, time = "year", spde = spde, family = tweedie(link = "log"),
                     epsilon_model = "loglinear")

sub$resid = sub$juv_cpue_kg_km2 - predict(juv_tweedie)$est

# juveniles


# adults


# save model
saveRDS(juv_tweedie, file=paste0(comm_name,"_nonstationary_juv.rds"))
