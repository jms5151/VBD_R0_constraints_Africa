# load libraries
library(rstan)

# load data (source = 'format_data_for_R0_stan_model.R)
load('../VBD-data/model_data_zikv.RData')

# avoid recompilation
rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

# fit model
stan_model_fit_zikv <- 
  sampling(
    stan_model('code/R0_model.stan')
  , data = model_data_zikv
  , iter = 10000
  , cores = parallel::detectCores()
  # , cores = 3
  )

# summary
stan_model_fit_zikv

# save stanfit object
saveRDS(stan_model_fit_zikv,'../models/stan_model_fit_zikv.rds')

