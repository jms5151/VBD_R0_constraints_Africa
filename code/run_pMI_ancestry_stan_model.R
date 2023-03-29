# load libraries
library(rstan)
library(matrixStats)
library(tidyverse)

# load data
zikv_afr_panel <- read.delim('../VBD-data/reformatted_ZIKV_african_panel.txt')

# Make population and virus factored data
zikv_afr_panel[, c('Population', 'Virus')] <- lapply(
  zikv_afr_panel[, c('Population', 'Virus')]
  , function(x) as.factor(x)
)

# Format data, wide to long
zikv_afr_panel <- zikv_afr_panel %>%
  gather(key = 'Infection', 'N', -c(Population, Dose, Virus, anc))

zikv_afr_panel <- zikv_afr_panel[rep(seq(nrow(zikv_afr_panel)), zikv_afr_panel$N),]

# make infection status 0/1
zikv_afr_panel$Infection <- ifelse(zikv_afr_panel$Infection == 'uninfected', 0, 1)

# for each trait, order by temperature from low to high 
zikv_afr_panel <- zikv_afr_panel %>% arrange(anc)

model_data <- 
  list(
    pMI_ancestry_N = nrow(zikv_afr_panel)
    , pMI_ancestry = zikv_afr_panel$anc
  )

# fit model
stan_model_fit_ancestry_pMI <- sampling(
  stan_model('code/model_infection_given_ancestry.stan')
  , data = model_data
  , iter = 1000
)

mod_params <- c('pMI_ancestry_intercept', 'beta')
rstan::traceplot(stan_model_fit_ancestry_pMI, par = c('lp__', mod_params), ncol = 2)
