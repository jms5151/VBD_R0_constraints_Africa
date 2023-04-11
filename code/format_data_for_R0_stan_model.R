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

# create simulated new data
temp_new  <- seq(10,40,0.1)

anc_new <- seq(0, 1, 0.1)
dose_new <- unname(quantile(zikv$Dose, c(0.25, 0.5, 0.75)))
virus_new <- c(0,1)

zikv_new <- expand.grid(anc_new, dose_new, virus_new)
colnames(zikv_new) <- c('anc', 'Dose', 'Virus')
zikv_new <- zikv_new %>% arrange(anc, Dose, Virus)

# survey data
survey_data <- read.csv('../VBD-data/combined_meta_allpops.csv')
survey_data <- subset(survey_data, !is.na(prop_aaa_ancestry))

survey_new <- expand.grid(survey_data$Location, dose_new, virus_new)
colnames(survey_new) <- c('Location', 'Dose', 'Virus')
survey_new <- survey_new %>%
  left_join(survey_data[,c('Location', 'prop_aaa_ancestry', 'bio.bio8_temp_wetq')])

# combine data for Zika model
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
    , climate_temp_new = temp_new                                                  # vector of temperatures to predict on
    , climate_N_new = length(temp_new)                                             # number of new temperatures
    , omega_ancestry_N = nrow(omega.df)
    , omega_ancestry = omega.df$prob
    , omega_ancestry_aa = omega.df$prop_aaa_ancestry
    , pMI_ancestry_N = nrow(zikv)
    , pMI_ancestry_Trials = zikv$Trials
    , pMI_ancestry_Infected = zikv$Infected
    , pMI_ancestry_propInf = zikv$Infected / zikv$Trials # only needed for assessing output
    , pMI_ancestry_X = model.matrix(~ 0 + scale(anc) + log(Dose), data = zikv)
    , pMI_ancestry_Virus = ifelse(zikv$Virus == 'ZIKV_Senegal_2011', 1, 0)
    , ancestry_N_new = nrow(zikv_new)
    , ancestry_X_new = model.matrix(~0 + scale(anc) + log(Dose), data = zikv_new)
    , ancestry_aa_new = zikv_new$anc
    , ancestry_Virus_new = zikv_new$Virus
    , ancestry_dose_new = zikv_new$Dose # only needed for assessing output
    , surveys_N = nrow(survey_new)
    , surveys_temp = survey_new$bio.bio8_temp_wetq
    , surveys_aa = survey_new$prop_aaa_ancestry
    , surveys_X = model.matrix(~0 + scale(prop_aaa_ancestry) + log(Dose), data = survey_new)
    , surveys_Virus = survey_new$Virus
    , surveys_dose = survey_new$Dose # only needed for assessing output
  )

save(model_data_zikv, file = '../VBD-data/model_data_zikv.RData')
