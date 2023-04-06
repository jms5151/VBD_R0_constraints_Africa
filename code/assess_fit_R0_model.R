# load libraries
library(rstan)
library(shinystan)
library(boot)
library(matrixStats)
library(tidyverse)
library(ggplot2)
library(cowplot)

# list parameters
# ignores sigmas, deltas, etc.
params <- c('omega_ancestry_constant'
            , 'omega_ancestry_d'
            , 'omega_ancestry_e'
            , 'alpha_climate_Tmin'
            , 'alpha_climate_Tmax'
            , 'alpha_climate_constant'
            , 'b_climate_Tmin'
            , 'b_climate_Tmax'
            , 'b_climate_constant'
            # , 'pMI_climate_Tmin'
            # , 'pMI_climate_Tmax'
            # , 'pMI_climate_constant'
            , 'pMI_climate_rmax'
            , 'pMI_climate_Topt'
            , 'pMI_climate_a'
            , 'pMI_ancestry_b0'
            , 'pMI_ancestry_beta[1]'
            , 'pMI_ancestry_beta[2]'
            , 'pMI_ancestry_gamma'
            , 'EIR_climate_Tmin'
            , 'EIR_climate_Tmax'
            , 'EIR_climate_constant'
            , 'lf_climate_Tmin'
            , 'lf_climate_Tmax'
            , 'lf_climate_constant' 
            )

# list prior data
prior_data <- c('alpha_climate'
                , 'b_climate'
                , 'pMI_climate'
                , 'EIR_climate'
                , 'lf_climate')

# open model
r0_mod <- readRDS('../models/stan_model_fit_zikv.rds')
load('../VBD-data/model_data_zikv.RData')
br_mod <- readRDS('../models/stan_model_fit_ancestry_alpha.rds')

# assess model fit
launch_shinystan(r0_mod)

# trace plots ------------------------------------------------------
pdf('figures/R0_stan_zikv_traceplots.pdf', width = 11, height = 8.5)
rstan::traceplot(r0_mod, par = c('lp__', params), ncol = 5)
dev.off()

# ppc plots --------------------------------------------------------

# extract ppc samples and summarise
ppc_indexes <- names(r0_mod)[grep('ppc', names(r0_mod))]
ppc_estimates <- rstan::extract(r0_mod, ppc_indexes)
ppc_estimates <- lapply(ppc_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)

ppc_estimates_quants <- do.call(rbind.data.frame, ppc_estimates)
colnames(ppc_estimates_quants) <- c('lower', 'median', 'upper')
ppc_estimates_quants$trait <- names(ppc_estimates)
# ppc_estimates_quants$trait <- gsub('_climate.*|_ppc.*|pMI_ancestry_', '', ppc_estimates_quants$trait)
ppc_estimates_quants$trait <- gsub('_ppc.*', '', ppc_estimates_quants$trait)

# load data
# mod_data <- readRDS('../VBD-Data/aa_traits_denv.RData')

prior_data_df <- mod_data[prior_data]
prior_data_df <- map_df(prior_data_df, ~as.data.frame(.x), .id='trait2')
colnames(prior_data_df)[2] <- 'value'

# combine estimated and observed data
ppc <- cbind(ppc_estimates_quants, prior_data_df)

# plot
pdf('figures/R0_ppc_plots.pdf', width = 11, height = 8.5)
ggplot(ppc, aes(value, median)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_classic() + 
  facet_wrap(.~trait, scales = 'free') +
  geom_abline()
dev.off()


# trait fits --------------------------------------------------
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
  plot(df[['climate_temp_new']], sampMeans, type='l', lwd=2, ylab=param_name, xlab=expression(paste("Temperature (",degree,"C)")), ylim = c(0,yMax), main = param_name)
  lines(df[['climate_temp_new']], sampQuantiles[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(df[['climate_temp_new']], sampQuantiles[,2], lty=2, col='red', ylim=c(0,yMax))
  if(param_name != 'R0'){
    param_name_orig <- paste0(param_name, '_climate_temp')
    param_vals <- paste0(param_name, '_climate')
    points(df[[param_name_orig]], df[[param_vals]], pch=16, ylim=c(0,yMax))
  }
}

pdf('figures/R0_Zika_trait_fit_plots.pdf', width = 11, height = 8.5)
par(mfrow = c(2, 3)) 
plotSamples(mod = r0_mod, param_name = 'alpha', df = mod_data)
plotSamples(mod = r0_mod, param_name = 'b', df = mod_data)
plotSamples(mod = r0_mod, param_name = 'pMI', df = mod_data)
plotSamples(mod = r0_mod, param_name = 'EIR', df = mod_data)
plotSamples(mod = r0_mod, param_name = 'lf', df = mod_data)
plotSamples(mod = r0_mod, param_name = 'R0', df = mod_data)
dev.off()

# prior vs posterior plots -----------------------------------

# function to overlay histogram and normal density
overlay_distributions_plot <- function(mod, param_name, priorMean, priorSD){
  # get posterior distribution
  df <- as.data.frame(rstan::extract(mod, param_name))
  # plot
  ggplot(df, aes_string(param_name)) +
    geom_histogram(aes(y = after_stat(density)), color = 'black', fill =  'yellow', alpha = 0.4) +
    stat_function(
      fun = dnorm, 
      args = list(mean = priorMean, sd = priorSD), 
      lwd = 2, 
      col = 'black'
    ) +
    theme_classic() +
    ggtitle(param_name) +
    ylab('') +
    xlab('')
}


a_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_constant', priorMean = 2.02E-04, priorSD = 0.01)
a_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmin', priorMean =  13.35, priorSD = 20)
a_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmax', priorMean = 40.08, priorSD = 20)
b_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_constant', priorMean = 8.49E-04, priorSD = 0.01)
b_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmin', priorMean = 17.05, priorSD = 20)
b_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmax', priorMean = 35.83, priorSD = 20)
pMI_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_constant', priorMean = 4.91E-04, priorSD = 0.01)
pMI_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Tmin', priorMean = 12.22, priorSD = 20)
pMI_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Tmax', priorMean = 37.46, priorSD = 20)
EIR_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_constant', priorMean = 6.65E-05, priorSD = 0.01)
EIR_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_Tmin', priorMean = 10.68, priorSD = 20)
EIR_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_Tmax', priorMean = 45.90, priorSD = 20)
lf_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_constant', priorMean = -1.48E-01, priorSD = 0.1)
lf_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_Tmin', priorMean = 9.16, priorSD = 20)
lf_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_Tmax', priorMean = 37.73, priorSD = 20)

p <- plot_grid(a_clim_const
          , a_clim_Tmin
          , a_clim_Tmax
          , b_clim_const
          , b_clim_Tmin
          , b_clim_Tmax
          , pMI_clim_const
          , pMI_clim_Tmin
          , pMI_clim_Tmax
          , EIR_clim_const
          , EIR_clim_Tmin
          , EIR_clim_Tmax
          , lf_clim_const
          , lf_clim_Tmin
          , lf_clim_Tmax
          )
ggsave('figures/prior_vs_posterior_plots.pdf', p, width = 11, height = 11)



# Compare parameter outputs with Mordecai PLoS NTD ----------------
list_of_draws <- rstan::extract(r0_mod)

# temperature at which R0 is optimized
r0Samples <- data.frame(list_of_draws['R0'])
r0Means<-colMeans(r0Samples, na.rm = T)
mod_data$climate_temp_new[which(r0Means == max(r0Means))]

param_names <- names(list_of_draws)[grepl('constant|min|max', names(list_of_draws))==T]
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

dodge <- position_dodge(width=0.2) 

pdf('figures/R0_trait_fit_compare_with_Mordecai_etal.pdf', width = 11, height = 6)
ggplot(param_vals2, aes(x = mean, y = parameter, color = study)) +
  geom_point(size = 2, position = dodge) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width=.05, position = dodge) +
  theme_bw() +
  facet_wrap(.~trait, ncol = 5) +
  coord_fixed(ratio = 60) +
  xlab('') +
  ylab('') +
  theme(legend.title = element_blank(), legend.position = 'top') +
  scale_colour_manual(values = c('darkblue','orange'))
dev.off()
