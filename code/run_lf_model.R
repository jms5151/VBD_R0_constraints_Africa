library(rstan)
library(parallel)


load('../VBD-data/model_data_zikv.RData')

new_list <- list('lf_climate_N' = model_data_zikv$lf_climate_N
                 , 'lf_climate_temp' = model_data_zikv$lf_climate_temp
                 , 'lf_climate' = model_data_zikv$lf_climate
                 , 'N_new' = model_data_zikv$N_new
                 , 'temp_new' = model_data_zikv$temp_new
) 

mod <- 
  sampling(
    stan_model('code/lf_model.stan')
    , data = new_list
    , iter = 1000
    , cores = parallel::detectCores()-1
    , chains = 3
  )

mod