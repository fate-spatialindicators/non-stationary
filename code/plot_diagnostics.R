library(dplyr)
library(sdmTMB)

species = read.csv("survey_data/species_list.csv", fileEncoding="UTF-8-BOM")
names(species) = tolower(names(species))
species = dplyr::rename(species,
                        common_name = common.name,
                        scientific_name = scientific.name)

.sp <- sub(" ", "_", species$common_name)
files <- paste0("output/", .sp, "_all_models.RData")
#files <- paste0("output/", .sp, "_ar1.RData")

for(i in 1:length(files)){

  load(files[i])
  ad_fit$resids <- residuals(ad_fit)
  ad_fit_ll$resids <- residuals(ad_fit_ll)

  jpeg(paste0("plots/diagnostics_",unique(ad_fit$data$common_name),".png"),
      width = 850, height = 500)
  par(mfrow = c(1, 2))

  q1 <- qqnorm(ad_fit$resids, plot.it = FALSE)
  q2 <- qqnorm(ad_fit_ll$resids, plot.it = FALSE)
  plot(range(q1$x, q2$x), range(q1$y, q2$y), type = "n",
       main = unique(ad_fit$data$common_name),
       xlab = "theoretical quantiles", ylab = "sample quantiles")
  abline(a = 0, b = 1)
  points(q1, col = "red")
  points(q2)
  legend("bottomright", legend = c("null", "loglinear"),
         pch = 1, col = c("red", "black"))

  h1 <- hist(ad_fit$resids, plot = FALSE)
  h2 <- hist(ad_fit_ll$resids, plot = FALSE)
  plot(h1, col=rgb(1,0,0,1/4), xlim = range(ad_fit$resids, ad_fit_ll$resids),
       main = NULL, xlab = "residuals")
  plot(h2, col=rgb(0,0,0,1/4), xlim = range(ad_fit$resids, ad_fit_ll$resids), add = TRUE)
  box()

  dev.off()
}

