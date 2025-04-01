
# LOAD DATA ------

library(tidyverse)

data <- read.csv(file = "otu_metadata_sy_20250318.csv")

binarise1 <- function(x, na.rm = FALSE) vegan::decostand(x, method = "pa")

# species to remove
spptoremove <- data %>%
  as.data.frame() %>%
  group_by(location.x) %>%
  summarise(across(c(17:412), function(x){
    sum(x > 0)
  })) %>%
  arrange(location.x) %>%
  dplyr::select(c(2:397)) %>%
  dplyr::mutate(across(c(1:396), binarise1)) |>
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) < 7)) %>%
  as.matrix %>%
  colnames()

OTUtable <- data %>%
  as.data.frame() %>%
  dplyr::filter(type == "L" | type == "W") %>%
  dplyr::select(c(17:412)) %>%
  dplyr::select(-one_of(spptoremove)) %>%
  as.matrix

data_infos <- data %>%
  dplyr::filter(type == "L" | type == "W") %>%
  rename(Site = location.x,
         Sample = sample,
         Primer = primer) %>%
  dplyr::select(Site, Sample, Primer, latitude, longitude, altitude, type, zone, NDVI,
                input_volume, elution_volume) %>%
  as.data.frame


data <- list(info = data_infos,
             OTU = OTUtable)

data$OTU[data$info$elution_volume == 0,] <- NA

# RUN MODEL -----

fitmodel  <- runOccPlus(data,
                        d = 2,
                        occCovariates = c(),
                        ordCovariates = c("altitude"),
                        detCovariates = c("input_volume"),
                        numSamples = 1000)

# OUTPUT  ---------

# plotOccupancyCovariates(fitmodel,
#                         covName = "")

plotDetectionCovariates(fitmodel,
                        covName = "input_volume")

plotOrdinationCovariates(fitmodel,
                         covName = "altitude")

plotOccupancyRates(fitmodel,
                   idx_species = 1:20)

plotCollectionRates(fitmodel,
                    idx_species = 1:20)

plotDetectionRates(fitmodel,
                   idx_species = 1:20)

plotFPDetectionRates(fitmodel,
                     idx_species = 1:5)

plotReadIntensity(fitmodel)
