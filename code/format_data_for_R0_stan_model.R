# format data
library(tidyverse)

# temperature dependent trait data ------------------------------------
traits.df <- read.csv('../VBD-Data/aegypti_traits_temp_formatted.csv')

# for each trait, order by temperature from low to high 
traits.df <- traits.df %>% arrange(trait_name_new, Temperature)

# ancestry dependent biting prob data ---------------------------------
omega.df <- read.csv('../VBD-data/combined_meta_allpops.csv')
omega.df <- subset(omega.df, !is.na(prop_aaa_ancestry))

# for each trait, order by temperature from low to high 
omega.df <- omega.df %>% arrange(prop_aaa_ancestry)

# ancestry dependent pMI data (combine 2 datasets) --------------------
zikv_afr_panel <- read.delim('../VBD-data/reformatted_ZIKV_african_panel.txt')

# format
zikv_afr_panel$Trials <- zikv_afr_panel$uninfected + zikv_afr_panel$infected
colnames(zikv_afr_panel)[5] <- 'Infected'

# add 'ZIKV' to virus name
zikv_afr_panel$Virus <- paste0('ZIKV_', zikv_afr_panel$Virus)

# load cape verde data
zikv_cpv <- read.csv('../VBD-data/CPV-ZIKV.csv')

zikv_cpv_wide <- zikv_cpv %>%
  group_by(Virus, Population, Titer) %>%
  summarise(Infected = sum(Infection == 1), Trials = length(Infection)) %>%
  as.data.frame()

# rename titer as dose
colnames(zikv_cpv_wide)[3] <- 'Dose'

# add ancestry
zikv_cpv_wide$anc <- 0.23 # CPV
zikv_cpv_wide$anc[zikv_cpv_wide$Population == 'NGO'] <- 0.3738158
zikv_cpv_wide$anc[zikv_cpv_wide$Population == 'Gabon'] <-  0.073
zikv_cpv_wide$anc[zikv_cpv_wide$Population == 'Guadeloupe'] <- 1

# combine
zikv_colnames <- c('Virus' ,'Population', 'anc', 'Dose', 'Infected', 'Trials')

zikv <- rbind(
  zikv_afr_panel[, zikv_colnames],
  zikv_cpv_wide[, zikv_colnames]
)

# order the data
zikv <- zikv %>% arrange(Trials, Infected)

# create simulated new datasets and format to combine --------------------------

# ancestry ---
anc_new <- seq(0, 1, 0.1)
dose_new <- unname(quantile(zikv$Dose, c(0.25, 0.5, 0.75)))
virus_new <- c(0,1)

zikv_new <- expand.grid(anc_new, dose_new, virus_new)
colnames(zikv_new) <- c('anc', 'Dose', 'Virus')

zikv_new <- zikv_new %>% 
  arrange(anc, Dose, Virus) %>%
  mutate(
    'Location' = rep(NA, nrow(zikv_new))
    , 'lon' = rep(NA, nrow(zikv_new))
    , 'lat' = rep(NA, nrow(zikv_new))
    , 'year' = rep(NA, nrow(zikv_new))
    , 'validation_type' = 'ancestry'
  )

zikv_new <- cbind('temp' = rep(29, nrow(zikv_new)), zikv_new)

# temperature ---
temp_new  <- seq(10,40,0.1)

temp_df <- data.frame(
  'temp' = temp_new
  , 'anc' = rep(0.5, length(temp_new))
  , 'Dose' = rep(min(zikv_new$Dose), length(temp_new))
  , 'Virus' = rep(0, length(temp_new))
  , 'Location' = rep(NA, length(temp_new))
  , 'lon' = rep(NA, length(temp_new))
  , 'lat' = rep(NA, length(temp_new))
  , 'year' = rep(NA, length(temp_new))
  , 'validation_type' = 'temperature'
)

# survey data ---
survey_data <- read.csv('../VBD-data/combined_meta_allpops.csv')
survey_data <- subset(survey_data, !is.na(prop_aaa_ancestry))

survey_new <- expand.grid(survey_data$Location, dose_new, virus_new)
colnames(survey_new) <- c('Location', 'Dose', 'Virus')
survey_new <- survey_new %>%
  left_join(survey_data[,c('Location', 'prop_aaa_ancestry', 'bio.bio8_temp_wetq')]) %>%
  mutate('lat' = rep(NA, nrow(survey_new))
         , 'lon' = rep(NA, nrow(survey_new))
         , 'year' = rep(NA, nrow(survey_new))
         , 'validation_type' = 'surveys'
         )

colnames(survey_new)[4:5] <- c('anc', 'temp')
survey_new <- survey_new[, colnames(temp_df)]

# map data ---
# map_data <- read.csv('../VBD-data/map_variables.csv')
# colnames(map_data) <- c('lon', 'lat', 'anc', 'temp')
# 
# map_data <- map_data %>%
#   mutate(
#     'Dose' = rep(min(zikv_new$Dose), nrow(map_data))
#     , 'Virus'= rep(0, nrow(map_data))
#     , 'validation_type' = 'map'
#     , 'Location' = NA
#   )
#   
# map_data <- map_data[, colnames(temp_df)]

# new_data <- new_data[1:1500, ]

# African cities time series data ---
# load & format aaa data
aaa_cities <- read.csv('../VBD-data/african_cities_2100aaa.csv')

# wide to long
aaa_cities <- aaa_cities %>% gather('Variable', 'anc', bio.bio15:tip) 

# remove variables other than IDs, ancestry, and year
aaa_cities <- aaa_cities[str_detect(aaa_cities$Variable, '^aaa'),]

# Make new Year column
aaa_cities$year <- as.numeric(gsub('aaa', '', aaa_cities$Variable))

# load temperature data
cmip6ssp585 <- read.csv('../VBD-data/CMIP6_TAS_Timeseries_GFDL-ESM4_ssp585.csv')

# format temperature from kelvin to celsius
cmip6ssp585$temp <- cmip6ssp585$mean - 273.15

# combine ancestry and temperature data for cities
big_cities <- aaa_cities %>% 
  left_join(cmip6ssp585[,c('City', 'year', 'temp')]) %>%
  mutate('Dose' = rep(min(zikv_new$Dose), nrow(big_cities))
    , 'Virus' = rep(0, nrow(big_cities))
    , 'Location' = paste0(City, ', ', Country)
    , 'lon' = Longitude
    , 'lat' = Latitude
    , 'validation_type' = 'big_cities'
  )

big_cities <- big_cities[, colnames(zikv_new)]

# combine new data ---
new_data <- do.call(rbind, list(
  zikv_new
  , temp_df
  , survey_new
  # , map_data
  , big_cities
))

# combine all data to fit and generate data with Zika model
model_data_zikv <-
  list(
    alpha_climate_N = sum(traits.df$trait_name_new == 'alpha')                           # number of observations
    , alpha_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'alpha']    # vector of temperatures
    , alpha_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'alpha']     # vector of trait values
    , b_climate_N = sum(traits.df$trait_name_new == 'b_zikv')                            # number of observations
    , b_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'b_zikv']       # vector of temperatures
    , b_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'b_zikv']        # vector of trait values
    , pMI_climate_N = sum(traits.df$trait_name_new == 'pMI_zikv')                        # number of observations
    , pMI_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'pMI_zikv']   # vector of temperatures
    , pMI_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'pMI_zikv']    # vector of trait values
    , EIR_climate_N = sum(traits.df$trait_name_new == 'EIR')                             # number of observations
    , EIR_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'EIR']        # vector of temperatures
    , EIR_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'EIR']         # vector of trait values
    , lf_climate_N = sum(traits.df$trait_name_new == 'lifespan')                         # number of observations
    , lf_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'lifespan']    # vector of temperatures
    , lf_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'lifespan']     # vector of trait values
    , climate_temp_new = temp_new                                                        # vector of temperatures to predict on
    , climate_N_new = length(temp_new)                                                   # number of new temperatures
    , omega_ancestry_N = nrow(omega.df)
    , omega_ancestry = omega.df$prob
    , omega_ancestry_aa = omega.df$prop_aaa_ancestry
    , pMI_ancestry_N = nrow(zikv)
    , pMI_ancestry_Trials = zikv$Trials
    , pMI_ancestry_Infected = zikv$Infected
    , pMI_ancestry_propInf = zikv$Infected / zikv$Trials # identifier for assessing output
    , pMI_ancestry_X = model.matrix(~ 0 + scale(anc) + log(Dose), data = zikv)
    , pMI_ancestry_Virus = ifelse(zikv$Virus == 'ZIKV_Senegal_2011', 1, 0)
    , N_new = nrow(new_data)
    , X_new = model.matrix(~0 + scale(anc) + log(Dose), data = new_data)
    , aa_new = new_data$anc
    , Virus_new = new_data$Virus
    , dose_new = new_data$Dose
    , temp_new = new_data$temp
    # identifiers for assessing output
    , location = new_data$Location
    , lon = new_data$lon
    , lat = new_data$lat
    , year = new_data$year
    , validationtype = new_data$validation_type
    # , ancestry_N_new = nrow(zikv_new)
    # , ancestry_X_new = model.matrix(~0 + scale(anc) + log(Dose), data = zikv_new)
    # , ancestry_aa_new = zikv_new$anc
    # , ancestry_Virus_new = zikv_new$Virus
    # , ancestry_dose_new = zikv_new$Dose # only needed for assessing output
    # , surveys_N = nrow(survey_new)
    # , surveys_temp = survey_new$bio.bio8_temp_wetq
    # , surveys_aa = survey_new$prop_aaa_ancestry
    # , surveys_X = model.matrix(~0 + scale(prop_aaa_ancestry) + log(Dose), data = survey_new)
    # , surveys_Virus = survey_new$Virus
    # , surveys_dose = survey_new$Dose # only needed for assessing output
    # , surveys_location = survey_new$Location # only needed for assessing output
    # , map_N = nrow(map_data)
    # , map_temp = map_data$wc2.1_10m_bio_8
    # , map_aa = map_data$prop_aaa_raster
    # , map_Virus = map_data$Virus
    # , map_X = model.matrix(~0 + scale(prop_aaa_raster) + log(Dose), data = map_data)
    # , map_lon = map_data$x # only needed for mapping output
    # , map_lat = map_data$y # only needed for mapping output
  )

save(model_data_zikv, file = '../VBD-data/model_data_zikv.RData')

