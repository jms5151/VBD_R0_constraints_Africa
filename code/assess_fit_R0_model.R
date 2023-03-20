# load libraries
library(rstan)
library(shinystan)
library(boot)
library(matrixStats)
library(tidyverse)

# list parameters
params <- c('alpha_climate_Tmin'
            , 'alpha_climate_Tmax'
            , 'alpha_climate_constant'
            , 'b_climate_Tmin'
            , 'b_climate_Tmax'
            , 'b_climate_constant'
            , 'pMI_climate_Tmin'
            , 'pMI_climate_Tmax'
            , 'pMI_climate_constant'
            , 'EIR_climate_Tmin'
            , 'EIR_climate_Tmax'
            , 'EIR_climate_constant'
            , 'lf_climate_Tmin'
            , 'lf_climate_Tmax'
            , 'lf_climate_constant')

# list prior data
prior_data <- c('alpha_climate'
                , 'b_climate'
                , 'pMI_climate'
                , 'EIR_climate'
                , 'lf_climate')

# open model
# r0_mod <- readRDS('../models/stan_model_fit_denv.rds')
r0_mod <- readRDS('../models/R0_stan_model_fit.rds')


# assess model fit
launch_shinystan(r0_mod)

# trace plots
pdf('figures/R0_stan_traceplots.pdf', width = 11, height = 8.5)
rstan::traceplot(r0_mod, par = c('lp__', params), ncol = 4)
dev.off()

# ppc

# extract ppc samples and summarise
ppc_indexes <- names(r0_mod)[grep('ppc', names(r0_mod))]
ppc_estimates <- rstan::extract(r0_mod, ppc_indexes)
ppc_estimates <- lapply(ppc_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)

ppc_estimates_quants <- do.call(rbind.data.frame, ppc_estimates)
colnames(ppc_estimates_quants) <- c('lower', 'median', 'upper')
ppc_estimates_quants$trait <- names(ppc_estimates)
ppc_estimates_quants$trait <- gsub('_climate.*', '', ppc_estimates_quants$trait)

# load data
traits.df <- read.csv('../VBD-Data/aegypti_traits_temp_formatted.csv')
traits.df$trait_name_new <- gsub('lifespan', 'lf', traits.df$trait_name_new)

# remove traits for disease not modeled
traits <- traits.df[!grepl('.*_zikv', traits.df$trait_name_new),]
traits$trait_name_new <- gsub('_denv', '', traits$trait_name_new) 

# for each trait, order by temperature from low to high 
traits <- traits %>% arrange(trait_name_new, Temperature)


# still not quite right -- not matching up,  temp out of order?
x <- cbind(ppc_estimates_quants, traits)
# library(purrr)
# 
mod_data <- readRDS('../VBD-Data/aa_traits_denv.RData')

prior_data_df <- mod_data[prior_data]
prior_data_df <- map_df(prior_data_df, ~as.data.frame(.x), .id='trait2')
colnames(prior_data_df)[2] <- 'value'
# 
ppc <- cbind(ppc_estimates_quants, prior_data_df)
# 
# library(ggplot2)
ggplot(x, aes(trait_value_new, median)) +
  geom_point() +
  theme_classic() +
  facet_wrap(.~trait, scales = 'free')
# 
# 
# test <- data.frame('trait' = mod_data$b_climate, 'temp' = mod_data$b_climate_temp)
# 
# test <- traits.df[traits.df$trait_name_new == 'b_denv',]
# 
# test2 <- cbind(ppc_estimates_quants, test)
# plot(test2$median, test2$trait_value_new, pch = 16, xlim = c(0,1), ylim = c(0,1))
