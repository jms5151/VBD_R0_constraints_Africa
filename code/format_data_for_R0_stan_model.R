# format data
library(tidyverse)
library(gdata)

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
temp_new  <- seq(10,40,0.5)

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

# African cities time series data ---
# load & format aaa data
# aaa_cities <- read.csv('../VBD-data/african_cities_1970_2100.csv')
# aaa_cities2 <- read.csv('../VBD-data/african_cities_1970_2015.csv')

aaa_cities <- read.csv('../VBD-data/african_cities_1970_SSP585_1970_2020_and_2090-2100.csv')
cityNames <- read.csv('../VBD-data/african_cities_latlon.csv')

# format names
aaa_cities$City <- cityNames$City
aaa_cities$Country <- cityNames$Country

# get multimodel mean aaa for 2100
aaa_cities$aaa2100 <- rowMeans(aaa_cities[,c('aaa2100_INM.CM4.8'
                                             , 'aaa2100_INM.CM4.8'
                                             , 'aaa2100_MIROC6'
                                             , 'aaa2100_MPI.ESM1.2.HR'
                                             , 'aaa2100_MRI.ESM2.0')])

# combine
# aaa_cities$aaa2015 <- aaa_cities2$aaa2015

# remove sites from Rose et al. 2015
# aaa_cities <- aaa_cities[29:nrow(aaa_cities),]

# wide to long
# aaa_cities <- aaa_cities %>% gather('Variable', 'anc', bio.bio15:aaa2015) 

aaa_cities <- aaa_cities[, c('City', 'Country', 'Longitude', 'Latitude', 'aaa1970', 'aaa2020', 'aaa2100')]
aaa_cities <- aaa_cities %>% gather('Variable', 'anc', aaa1970:aaa2100) 

# remove variables other than IDs, ancestry, and year
aaa_cities <- aaa_cities[str_detect(aaa_cities$Variable, '^aaa'),]

# Make new Year column
aaa_cities$year <- as.numeric(gsub('aaa', '', aaa_cities$Variable))
# change end of century year to group
aaa_cities$year <- ifelse(aaa_cities$year > 2080, '2090-2100', aaa_cities$year) 

# average ancestry for end of century estimate
aaa_cities <- aaa_cities %>%
  filter(year == '1970' | year == '2020' | year == '2090-2100') %>%
  group_by(City, Country, Longitude, Latitude, year) %>%
  summarise(anc = mean(anc)) %>%
  as.data.frame()

# load temperature data
gfdl <- read.csv('../VBD-data/CMIP6_TAS_Timeseries_GFDL-ESM4_1970_2020_2090-2100.csv')
miroc <- read.csv('../VBD-data/CMIP6_TAS_Timeseries_MIROC6_1970_2020_2090-2100.csv')
NorESM2 <- read.csv('../VBD-data/CMIP6_TAS_Timeseries_NorESM2-MM_1970_2020_2090-2100.csv')
GISS <- read.csv('../VBD-data/CMIP6_TAS_Timeseries_GISS-E2-1-G_1970_2020_2090-2100.csv')
INM <- read.csv('../VBD-data/CMIP6_TAS_Timeseries_INM-CM4-8_1970_2020_2090-2100.csv')

# combine
cmip <- combine(gfdl, miroc, NorESM2, GISS, INM)

# format temperature from kelvin to celsius
cmip$temp <- cmip$mean - 273.15

# change end of century year to group
cmip$year <- ifelse(cmip$year > 2080, '2090-2100', cmip$year)

# average temperature year grouping
cmip_mean <- cmip %>%
  group_by(City, year) %>%
  summarise(temp = mean(temp)) %>%
  as.data.frame()

# fix city spelling to match aaa
# this is a terrible way of doing it, but it works
cmip_mean$City <- gsub('.\\¿½', 'é', cmip_mean$City)
cmip_mean$City <- gsub('_', ' ', cmip_mean$City)
cmip_mean$City[cmip_mean$City == 'Durban(eThekwini)'] <- sort(setdiff(aaa_cities$City, cmip_mean$City))[1]
cmip_mean$City[cmip_mean$City == 'East Rand(Ekurhuleni)'] <- sort(setdiff(aaa_cities$City, cmip_mean$City))[1]
cmip_mean$City[cmip_mean$City == 'Pretoria(Tshwane)'] <- sort(setdiff(aaa_cities$City, cmip_mean$City))[1]

# combine ancestry and temperature data for cities
big_cities <- aaa_cities %>% 
  left_join(cmip_mean[,c('City', 'year', 'temp')]) %>%
  mutate('Dose' = rep(min(zikv_new$Dose), nrow(aaa_cities))
    , 'Virus' = rep(0, nrow(aaa_cities))
    , 'Location' = paste0(City, ', ', Country)
    , 'lon' = Longitude
    , 'lat' = Latitude
    , 'validation_type' = 'big_cities'
  )

big_cities <- big_cities[, colnames(zikv_new)]

# contour data
temps <- seq(from = 15, to = 40, by =  1)
aaa <- seq(from = 0, to = 1, by = 0.1)

contour <- expand.grid('temp' = temps, 'anc' = aaa)

contour <- contour %>%
  mutate(
    'Dose' = rep(min(zikv_new$Dose), nrow(contour))
    , 'Virus' = rep(0, nrow(contour))
    , 'Location' = rep(NA, nrow(contour))
    , 'lon' = rep(NA, nrow(contour))
    , 'lat' = rep(NA, nrow(contour))
    , 'year' = rep(NA, nrow(contour))
    , 'validation_type' = 'contour'
  )

# combine new data ---
new_data <- do.call(rbind, list(
  zikv_new
  , temp_df
  , survey_new
  , big_cities
  , contour
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
  )

save(model_data_zikv, file = '../VBD-data/model_data_zikv.RData')

