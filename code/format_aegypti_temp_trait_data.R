# Format temperature experiments for Aedes aegypti and denv/zikv

# format data from Mordecai et al. 2017 PLoS NTD
# read data
aegypti_traits <- read.csv('../VBD-Data/trait_data/aegyptiDENVmodelTempData_2016-03-30.csv') 
aegypti_prior_traits <- read.csv('../VBD-Data/trait_data/Aedes_prior_data.csv') 

# rename temperature column
aegypti_traits$Temperature <- aegypti_traits$T

# format trait names and values
# we need a, b, c (pMI), mu, PDR (EIR)
aegypti_traits$trait_name_new <- aegypti_traits$trait.name
aegypti_traits$trait_value_new <- aegypti_traits$trait

# alpha (biting rate)
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == 'GCR'] <- 'alpha'

# b (probability mosquito infectious)
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == 'b' &
                                aegypti_traits$ref != 'Lambrects_et_al_2011_PNAS'] <- 'b_denv'

# lifespan (1/lifespan = mortality rate, mu)
aegypti_traits$trait_value_new[aegypti_traits$trait.name == 'p/days'] <- 1/aegypti_traits$trait[aegypti_traits$trait.name == 'p/days'] 
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == 'p/days'] <- 'lifespan'
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == '1/mu'] <- 'lifespan'

# pMI (probability of mosquito infectiousness)
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == 'c' &
                                aegypti_traits$ref != 'Lambrects_et_al_2011_PNAS'] <- 'pMI_denv'

# EIR
aegypti_traits$trait_value_new[aegypti_traits$trait.name == 'EIP'] <- 1/aegypti_traits$trait[aegypti_traits$trait.name == 'EIP'] 
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == 'EIP'] <- 'EIR'
aegypti_traits$trait_name_new[aegypti_traits$trait_name_new == 'PDR' &
                                aegypti_traits$ref != 'Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper'] <- 'EIR'

# keep only data needed for next gen R0 model
aa_traits_denv <- aegypti_traits[, c('Temperature', 'trait_name_new', 'trait_value_new', 'ref')]
aa_traits_denv <- subset(aa_traits_denv, trait_name_new == 'alpha'|
                      trait_name_new == 'lifespan'|
                      trait_name_new == 'b_denv'|
                      trait_name_new == 'pMI_denv'|
                      trait_name_new == 'EIR') 

## didn't seem to help
# # add prior data for aedes aegypti from Mordecai S2
# aa_prior <- subset(aegypti_prior_traits, trait.name == 'EIP')
# aa_prior$Temperature <- aa_prior$T
# aa_prior$trait_name_new <- 'EIR'
# aa_prior$trait_value_new <- 1/aa_prior$trait  
# aa_prior <- aa_prior[, c('Temperature', 'trait_name_new', 'trait_value_new', 'ref')]
# 
# aa_traits_denv <- rbind(aa_traits_denv, aa_prior)

# format data from Tesla et al. Proc B 2018
# load library
library(tidyverse)

# load data
aegypti_traits_zikv <- read.csv('../Tesla_2018_ProcB_data/InfectionData.csv') 

# format data
aegypti_traits_zikv$Temperature <- as.numeric(aegypti_traits_zikv$Temperature)

aa_traits_zikv <- aegypti_traits_zikv %>%
  group_by(Temperature, Batch) %>%
  summarise(
    b_zikv = sum(Infected == 1)/length(Batch)
    , pMI_zikv = sum(Infectious == 1)/length(Batch)
  ) %>%
  select(-Batch)  %>%
  gather('trait_name_new', 'trait_value_new', -Temperature) %>%
  mutate(ref = 'Tesla_2018_ProcB') %>%
  as.data.frame()


#  combine denv and zikv data
aa_traits <- rbind(aa_traits_denv, aa_traits_zikv)

# save data
write.csv(aa_traits, '../VBD-Data/aegypti_traits_temp_formatted.csv', row.names = F)
