# load libraries
library(rstan)
library(shinystan)
library(boot)
library(matrixStats)
library(tidyverse)

# load data
traits.df <- read.csv('../VBD-Data/aegypti_traits_temp_formatted.csv')

# for each trait, order by temperature from low to high 
traits.df <- traits.df %>% arrange(trait_name_new, Temperature)

# Dengue
# format data
model_data_denv <-
  list(
    alpha_climate_N = sum(traits.df$trait_name_new == 'alpha'),                         # number of observations
    alpha_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'alpha'],    # vector of temperatures
    alpha_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'alpha'],     # vector of trait values
    b_climate_N = sum(traits.df$trait_name_new == 'b_denv'),                            # number of observations
    b_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'b_denv'],       # vector of temperatures
    b_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'b_denv'],        # vector of trait values
    pMI_climate_N = sum(traits.df$trait_name_new == 'pMI_denv'),                        # number of observations
    pMI_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'pMI_denv'],   # vector of temperatures
    pMI_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'pMI_denv'],    # vector of trait values
    EIR_climate_N = sum(traits.df$trait_name_new == 'EIR'),                             # number of observations
    EIR_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'EIR'],        # vector of temperatures
    EIR_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'EIR'],         # vector of trait values
    lf_climate_N = sum(traits.df$trait_name_new == 'lifespan'),                         # number of observations
    lf_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'lifespan'],    # vector of temperatures
    lf_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'lifespan'],     # vector of trait values
    climate_temp_new = seq(10,40,0.1),                                                  # vector of temperatures to predict on
    climate_N_new = length(seq(10,40,0.1))                                              # number of new temperatures
  )

# save data
saveRDS(model_data_denv, '../VBD-Data/aa_traits_denv.RData')

# fit model
stan_model_fit_denv <- sampling(
  stan_model('code/R0_model.stan')
  , data = model_data_denv
  , iter = 10000
)

# save and open stanfit object
saveRDS(stan_model_fit_denv,'../models/stan_model_fit_denv.rds')
stan_model_fit_denv <- readRDS('../models/stan_model_fit_denv.rds')

# Zika Model isn't fitting
# If we need separate Zika and dengue models, suggest making sd of priors tighter (maybe 10 instead of 20)
# format data
# try fitting one zika variable at a time, b and pMI seem to have different issues
model_data_zikv <-
  list(
    alpha_climate_N = sum(traits.df$trait_name_new == 'alpha'),                         # number of observations
    alpha_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'alpha'],    # vector of temperatures
    alpha_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'alpha'],     # vector of trait values
    b_climate_N = sum(traits.df$trait_name_new == 'b_zikv'),                            # number of observations
    b_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'b_zikv'],       # vector of temperatures
    b_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'b_zikv'],        # vector of trait values
    pMI_climate_N = sum(traits.df$trait_name_new == 'pMI_zikv'),                        # number of observations
    pMI_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'pMI_zikv'],   # vector of temperatures
    pMI_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'pMI_zikv'],    # vector of trait values
    EIR_climate_N = sum(traits.df$trait_name_new == 'EIR'),                             # number of observations
    EIR_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'EIR'],        # vector of temperatures
    EIR_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'EIR'],         # vector of trait values
    lf_climate_N = sum(traits.df$trait_name_new == 'lifespan'),                         # number of observations
    lf_climate_temp = traits.df$Temperature[traits.df$trait_name_new == 'lifespan'],    # vector of temperatures
    lf_climate = traits.df$trait_value_new[traits.df$trait_name_new == 'lifespan'],     # vector of trait values
    climate_temp_new = seq(10,40,0.1),                                                  # vector of temperatures to predict on
    climate_N_new = length(seq(10,40,0.1))                                              # number of new temperatures
  )

# fit model
stan_model_fit_zikv <- sampling(
  # stan_model('code/zikv_stan_test_model.stan')
  stan_model('code/R0_model.stan')
  , data = model_data_zikv
  , iter = 8000
)


# save and open stanfit object
saveRDS(stan_model_fit_zikv,'../models/stan_model_fit_zikv.rds')
stan_model_fit_zikv <- readRDS('../models/stan_model_fit_zikv.rds')

plotSamples(mod = stan_model_fit_zikv, param_name = 'pMI', df = model_data_zikv)
plotSamples(mod = stan_model_fit_zikv, param_name = 'b', df = model_data_zikv)

# tpc <- function(c, Tmin, Tmax, m){
#   t_seq <- seq(10,45,0.1)
#   t_vec <- c()
#   for(i in t_seq){
#     t_vec <- c(t_vec, c*i*(i-Tmin)*sqrt(Tmax-i))
#   }
#   plot(t_seq, t_vec, type = 'l', ylim = c(0,0.30), main = paste0('c = ', c, ', Tmin = ', Tmin, ', Tmax = ', Tmax))
#   points(model_data_zikv$pMI_climate_temp, model_data_zikv$pMI_climate, pch = 16)
# }
# 
# tpc(c= 0.0002, 16, 38)

# try https://padpadpadpad.github.io/rTPC/articles/rTPC.html

# plots
rstan::traceplot(stan_model_fit_zikv, par = c('lp__', 'alpha_climate_Tmin','alpha_climate_Tmax','alpha_climate_constant'), ncol = 2)
rstan::traceplot(stan_model_fit_zikv, par = c('lp__', 'b_climate_Tmin','b_climate_Tmax','b_climate_constant'), ncol = 2)
rstan::traceplot(stan_model_fit_zikv, par = c('lp__', 'pMI_climate_Tmin','pMI_climate_Tmax','pMI_climate_constant'), ncol = 2)
rstan::traceplot(stan_model_fit_zikv, par = c('lp__', 'EIR_climate_Tmin','EIR_climate_Tmax','EIR_climate_constant'), ncol = 2)
rstan::traceplot(stan_model_fit_zikv, par = c('lp__', 'lf_climate_Tmin','lf_climate_Tmax','lf_climate_constant'), ncol = 2)
# bayesplot::mcmc_acf(as.matrix(stan_model_fit_denv), pars = c("b_climate_Tmin","b_climate_Tmax","b_climate_constant"))
# bayesplot::mcmc_areas(as.matrix(stan_model_fit_denv), pars = c("b_climate_Tmin","b_climate_Tmax","b_climate_constant"), prob = 0.95)

plotSamples <- function(mod, param_name, df){
  list_of_draws <- rstan::extract(mod)
  if(param_name != 'R0'){
    param_name_new <- paste0(param_name, '_climate_new')
  } else {
    param_name_new <- param_name
  }
  samps <- data.frame(list_of_draws[param_name_new])
  sampMeans <- colMeans(samps, na.rm = T)
  sampQuantiles <- colQuantiles(list_of_draws[[param_name_new]], na.rm = T, probs = c(0.025, 0.975))
  sampQuantiles <- ifelse(sampQuantiles < 0, 0, sampQuantiles)
  
  yMin = 0
  yMax = max(sampQuantiles)
  plot(df[['climate_temp_new']], sampMeans, type='l', lwd=2, ylab=param_name, xlab=expression(paste("Temperature (",degree,"C)")), ylim = c(0,yMax))
  lines(df[['climate_temp_new']], sampQuantiles[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(df[['climate_temp_new']], sampQuantiles[,2], lty=2, col='red', ylim=c(0,yMax))
  if(param_name != 'R0'){
    param_name_orig <- paste0(param_name, '_climate_temp')
    param_vals <- paste0(param_name, '_climate')
    points(df[[param_name_orig]], df[[param_vals]], pch=16, ylim=c(0,yMax))
  }
}

plotSamples(mod = stan_model_fit_denv, param_name = 'alpha', df = model_data_denv)
plotSamples(mod = stan_model_fit_denv, param_name = 'b', df = model_data_denv)
plotSamples(mod = stan_model_fit_denv, param_name = 'pMI', df = model_data_denv)
plotSamples(mod = stan_model_fit_denv, param_name = 'EIR', df = model_data_denv)
plotSamples(mod = stan_model_fit_denv, param_name = 'lf', df = model_data_denv)
plotSamples(mod = stan_model_fit_denv, param_name = 'R0', df = model_data_denv)

# parameter outputs, compare with Mordecai PLoS NTD
list_of_draws <- rstan::extract(stan_model_fit_denv)

# temperature at which R0 is optimized
r0Samples <- data.frame(list_of_draws['R0'])
r0Means<-colMeans(r0Samples, na.rm = T)
model_data_denv$climate_temp_new[which(r0Means == max(r0Means))]

param_names <- names(list_of_draws)[grepl('constant|min|max|opt|_a$', names(list_of_draws))==T]
param_means <- as.data.frame(lapply(list_of_draws[param_names], mean, na.rm = T))
param_quant <- lapply(list_of_draws[param_names], quantile, c(0.025, 0.975), na.rm = T)
param_quant <- as.data.frame(param_quant)

param_vals <- rbind(param_means, param_quant) 
param_vals$type <- row.names(param_vals)
param_vals <- param_vals %>% 
  gather('parameter', 'value', -type) %>%
  separate(parameter, c('trait', 'toremove', 'parameter')) %>%
  spread(type, value) %>%
  select(-toremove) %>%
  mutate(study = 'This study')
colnames(param_vals)[3:5] <- c('mean', 'lower', 'upper')

# combine with Mordecai
mordecai <- read.csv('../VBD-Data/trait_data/table_b_trait_values.csv')

param_vals2 <- rbind(param_vals, mordecai)
param_vals2[,c('trait', 'parameter')] <- lapply(param_vals2[,c('trait', 'parameter')], as.factor)

library(ggplot2)

dodge <- position_dodge(width=0.2) 

ggplot(param_vals2, aes(x = mean, y = parameter, color = study)) +
  geom_point(size = 2, position = dodge) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width=.05, position = dodge) +
  theme_bw() +
  # facet_wrap(.~trait+parameter, scales = 'free', ncol = 5) +
  facet_wrap(.~trait, ncol = 5) +
  # facet_grid(.~trait+parameter, cols = 5) +
  coord_fixed(ratio = 60) +
  xlab('') +
  ylab('') +
  theme(legend.title = element_blank(), legend.position = 'top') +
  scale_colour_manual(values = c('darkblue','orange'))

## other figures of interest
# Comparing ppc - data vs predictions
# comparing prior distributions vs posterior samples for all parameters