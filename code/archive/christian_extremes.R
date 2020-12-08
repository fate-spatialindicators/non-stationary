# this is all derived from code/run_wc_models.r (in fate-spatialindicators/consonants github repo)
#devtools::install_github("pbs-assess/sdmTMB")
library(sdmTMB)
library(dplyr)
library(sp)
library(nlme)

dat = readRDS("C:/Users/cjcco/OneDrive/Desktop/extremes/joined_nwfsc_data.rds")

dplyr::group_by(dat, species) %>% 
  summarize(min = min(latitude_dd[which(cpue_kg_km2 > 0)]),
            max = max(latitude_dd[which(cpue_kg_km2 > 0)])) %>% 
  as.data.frame() %>% arrange(min)

# UTM transformation
dat_ll = dat
coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
# convert to utm with spTransform
dat_utm = spTransform(dat_ll, 
                      CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
# convert back from sp object to data frame
dat = as.data.frame(dat_utm)
dat = dplyr::rename(dat, longitude = longitude_dd, 
                    latitude = latitude_dd)

df = expand.grid("species" = unique(dat$species),
                 spatial_only=c(FALSE), 
                 depth_effect = c(TRUE),
                 time_varying = c(FALSE)
)

for(i in 1:nrow(df)) {
  
  # filter by species, and select range within occurrences
  sub = dplyr::filter(dat, 
                      species == df$species[i])# %>% 
  #dplyr::filter(latitude > min(latitude[which(cpue_kg_km2>0)]),
  #  latitude <= max(latitude[which(cpue_kg_km2>0)]),
  #  longitude > min(longitude[which(cpue_kg_km2>0)]),
  #  longitude < max(longitude[which(cpue_kg_km2>0)]))
  
  # rescale variables
  sub$depth = scale(log(sub$depth))
  sub$o2 = scale(log(sub$o2))
  sub$temp = scale(sub$temp)
  
  # drop points with missing values
  sub = dplyr::filter(sub,!is.na(temp),!is.na(depth))
  
  # rename variables to make code generic
  sub = dplyr::rename(sub, enviro = as.character(df$covariate[i]))
  
  # make spde
  spde <- try(make_spde(x = sub$longitude, y = sub$latitude, 
                        n_knots = 150), silent=TRUE)
  if(class(spde) != "try-error") {
    formula = paste0("cpue_kg_km2 ~ -1")
    
    time_formula = "~ -1"
    
    #formula = paste0(formula, " + ", 
    #  "enviro", " + I(","enviro","^2)")
    time_varying = NULL
    time = NULL
    
    formula = paste0(formula, " + as.factor(year)")
    
    if(df$depth_effect[i]==TRUE) {
      formula = paste0(formula, " + depth + I(depth^2)")
    }
    
    # fit model
    m <- try(sdmTMB(
      formula = as.formula(formula),
      time_varying = NULL,
      spde = spde,
      time = "year",
      family = tweedie(link = "log"),
      data = sub,
      anisotropy = TRUE,
      spatial_only = FALSE
    ), silent=TRUE)
    
    #sd_report <- summary(m$sd_report)
    #params <- as.data.frame(sd_report[grep("quadratic", row.names(sd_report)), ])
    
    if(class(m)!="try-error") {
      
      # pull out epsilon_st, which is the spatiotemporal variation
      epsilon_st = summary(m$sd_report)[grep("epsilon_st",rownames(summary(m$sd_report))),]
      epsilon_st = as.data.frame(epsilon_st)
      locs = spde$mesh$loc
      epsilon_st$time = sort(rep(seq(1,nrow(epsilon_st)/nrow(locs)),
                                 nrow(locs)))
      epsilon_st$loc = rep(seq(1,nrow(locs)), length(unique(epsilon_st$time)))
      epsilon_st$year = unique(dat$year)[epsilon_st$time]
      epsilon_st$lon = locs[epsilon_st$loc,1]
      epsilon_st$lat = locs[epsilon_st$loc,2]      
      
      # 1. calculate variance for each year
      var = group_by(epsilon_st,year) %>% 
        summarize(v = var(Estimate))
      # 2. calculate variance on first-differenced fields
      var_diff = epsilon_st %>% 
        group_by(loc) %>% 
        dplyr::mutate(diff_est = c(NA,diff(Estimate))) %>% 
        ungroup() %>% 
        group_by(year) %>% 
        summarize(v_diff = var(diff_est,na.rm=T))
      # 3. calculate spatial autocorrelation by year
      var = dplyr::left_join(var,var_diff)
      
      var$range = NA
      var$nugget = NA
      for(j in 1:nrow(var)) {
        fit <- nlme::gls(Estimate ~ 1, 
                         data = dplyr::filter(epsilon_st,year==unique(epsilon_st$year)[j]),
                         correlation = corGaus(form = ~lat + lon, nugget = TRUE),
                         method = "REML")
        var$range[j] = attr(fit$apVar,"Pars")[1]
        var$nugget[j] = attr(fit$apVar,"Pars")[2]
      }
      
      # 3. calculate spatial autocorrelation by year
      var = dplyr::left_join(var,var_diff)
      
      var$diff_range = NA
      var$diff_nugget = NA
      var_diff = epsilon_st %>% 
        group_by(loc) %>% 
        dplyr::mutate(diff_est = c(NA,diff(Estimate)))
      for(j in 2:nrow(var)) {
        fit <- nlme::gls(diff_est ~ 1, 
                         data = dplyr::filter(var_diff,year==unique(epsilon_st$year)[j]),
                         correlation = corGaus(form = ~lat + lon, nugget = TRUE),
                         method = "REML")
        if(!is.null(attr(fit$apVar,"Pars"))) {
          var$diff_range[j] = attr(fit$apVar,"Pars")[1]
        var$diff_nugget[j] = attr(fit$apVar,"Pars")[2]
        }
      }
      saveRDS(var,file=paste0("C:/Users/cjcco/OneDrive/Desktop/extremes/summary_extremes_",i,".rds"))
      saveRDS(m, file=paste0("C:/Users/cjcco/OneDrive/Desktop/extremes/model_extremes_",i,".rds"))
      #sd_report <- summary(m$sd_report)
      #params <- as.data.frame(sd_report[grep("quadratic", row.names(sd_report)), ])
    }
  } # end try on spde
  
}
