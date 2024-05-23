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

# ancestry dependent pMI based on OID50 -------------------------------
oid50.df <- read.csv('../VBD-data/OID50s_ZIKV_Cambodia.csv')

# calculate pMI with reference to Guadeloupe strain
oid50.df$pMI <- ifelse(oid50.df$Study == 'Aubry et al.', oid50.df$OID50[oid50.df$Population == 'Guadeloupe' & oid50.df$Study == 'Aubry et al.'] / oid50.df$OID50, 
                       oid50.df$OID50[oid50.df$Population == 'Guadeloupe' & oid50.df$Study == 'Rose et al.'] / oid50.df$OID50)

# assume probability of most infectious strain is 50%, reduce all pMI values by 0.5 
oid50.df$pMI <- oid50.df$pMI - 0.70

# transform from 0-100 to 0-1 to match omega
oid50.df$Aaa <- oid50.df$Aaa/100

# replace 0% ancestry with extremely low value (model evaluates to inf at exactly 0)
oid50.df$Aaa[oid50.df$Aaa == 0] <- 0.000001

# order by ancestry
oid50.df <- oid50.df %>% arrange(Aaa)

# library(ggplot2)
# ggplot(oid50.df, aes(x = Aaa, y = pMI)) + geom_point(size = 6, alpha = 0.3) + theme_bw() + geom_smooth(method = 'lm', formula = y ~ poly(x, 2))

# create simulated new datasets and format to combine --------------------------

# ancestry ---
anc_new <- data.frame('anc' = seq(0, 1, 0.1))

anc_df <- anc_new %>% 
  mutate(
    'temp' = rep(29, nrow(anc_new))
    , 'Location' = rep(NA, nrow(anc_new))
    , 'lon' = rep(NA, nrow(anc_new))
    , 'lat' = rep(NA, nrow(anc_new))
    , 'year' = rep(NA, nrow(anc_new))
    , 'validation_type' = 'ancestry'
  )


# temperature ---
temp_new  <- data.frame('anc' = 0.5, 'temp' = seq(10, 40, 0.5))

temp_df <- temp_new %>% 
  mutate('Location' = rep(NA, nrow(temp_new))
    , 'lon' = rep(NA, nrow(temp_new))
    , 'lat' = rep(NA, nrow(temp_new))
    , 'year' = rep(NA, nrow(temp_new))
    , 'validation_type' = 'temperature'
  )

# survey data ---
survey_data <- read.csv('../VBD-data/combined_meta_allpops.csv')
survey_data <- subset(survey_data, !is.na(prop_aaa_ancestry))

# rename columns
survey_data <- survey_data %>%
  rename('anc' = 'prop_aaa_ancestry'
         , 'temp' = 'bio.bio8_temp_wetq')

survey_new <- survey_data[, c('anc', 'temp', 'Location')] %>%
  mutate('lat' = rep(NA, nrow(survey_data))
         , 'lon' = rep(NA, nrow(survey_data))
         , 'year' = rep(NA, nrow(survey_data))
         , 'validation_type' = 'surveys'
         )

# African cities time series data ---
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

# wide to long
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

# #IQR for paper
# cmip_mean <- cmip %>%
#   group_by(source, City, year) %>%
#   summarise(temp = mean(temp)) %>%
#   group_by(City, year) %>%
#   summarise(mtemp = mean(temp), minT = min(temp), maxT = max(temp), diff = maxT - minT) %>%
#   as.data.frame()
# quantile(cmip_mean$diff)
# 
# quantile(cmip_mean$diff[cmip_mean$year == '2020'])
# quantile(cmip_mean$diff[cmip_mean$year == '2090-2100'])

# combine ancestry and temperature data for cities
big_cities <- aaa_cities %>% 
  left_join(cmip_mean[,c('City', 'year', 'temp')]) %>%
  mutate('Location' = paste0(City, ', ', Country)
    , 'lon' = Longitude
    , 'lat' = Latitude
    , 'validation_type' = 'big_cities'
  )

big_cities <- big_cities[, colnames(survey_new)]

# contour data
aaa <- seq(from = 0, to = 1, by = 0.1)
temps <- seq(from = 15, to = 40, by =  1)

contour <- expand.grid('anc' = aaa, 'temp' = temps)

contour <- contour %>%
  mutate('Location' = rep(NA, nrow(contour))
    , 'lon' = rep(NA, nrow(contour))
    , 'lat' = rep(NA, nrow(contour))
    , 'year' = rep(NA, nrow(contour))
    , 'validation_type' = 'contour'
  )

# combine new data ---
new_data <- do.call(rbind, list(
  anc_df
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
    , climate_temp_new = temp_new$temp                                                   # vector of temperatures to predict on
    , climate_N_new = length(temp_new)                                                   # number of new temperatures
    , omega_ancestry_N = nrow(omega.df)
    , omega_ancestry = omega.df$prob
    , omega_ancestry_aa = omega.df$prop_aaa_ancestry
    , pMI_ancestry_N = nrow(oid50.df)
    , pMI_ancestry = oid50.df$pMI
    , pMI_ancestry_aa = oid50.df$Aaa 
    , N_new = nrow(new_data)
    , aa_new = new_data$anc
    , temp_new = new_data$temp
    # identifiers for assessing output
    , location = new_data$Location
    , lon = new_data$lon
    , lat = new_data$lat
    , year = new_data$year
    , validationtype = new_data$validation_type
  )

save(model_data_zikv, file = '../VBD-data/model_data_zikv.RData')

