library(dplyr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
future::plan(future::multisession, workers = max(floor(future::availableCores() / 4), 4L))

species <- read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) <- tolower(names(species))
species <- rename(species, common_name = common.name)
files <- paste0("output/", sub(" ", "_", species$common_name), "_ar1.RData")
# files <- paste0("output/", sub(" ", "_", species$common_name), "_all_models.RData")
models <- furrr::future_map(files, function(f) {
  cat(f, "\n")
  load(f)
  ad_fit_ll
})

get_b_eps <- function(x) {
  est <- as.list(x$sd_report, "Estimate", report = TRUE)$b_epsilon_logit
  se <- as.list(x$sd_report, "Std. Error", report = TRUE)$b_epsilon_logit
  tibble(
    est_b_eps = est,
    lwr_b_eps = est - 1.96 * se,
    upr_b_eps = est + 1.96 * se
  )
}

names(models) <- species$common_name
b_eps <- furrr::future_map_dfr(models, get_b_eps, .id = "common_name")
future::plan(future::sequential)

g1 <- b_eps %>%
  # filter(est_b_eps > -0.3, est_b_eps < 0.3) %>%
  ggplot(aes(est_b_eps, y = forcats::fct_reorder(common_name, est_b_eps),
    xmin = lwr_b_eps, xmax = upr_b_eps)) +
  geom_pointrange() +
  geom_vline(xintercept = 0, lty = 2) +
  ylab("")
g1
ggsave("plots/ar1-eps-time-effect.png", width = 7, height = 7)
