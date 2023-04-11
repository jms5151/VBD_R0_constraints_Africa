# load libraries
library(rstan)
# library(shinystan)
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

# open model
r0_mod <- readRDS('../models/stan_model_fit_zikv.rds')
load('../VBD-data/model_data_zikv.RData')
mod_data <- model_data_zikv

# assess model fit
# launch_shinystan(r0_mod)

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
ppc_estimates_quants$trait <- gsub('_ppc.*', '', ppc_estimates_quants$trait)

# load data
# mod_data <- readRDS('../VBD-Data/aa_traits_denv.RData')

prior_data_df <- mod_data[unique(ppc_estimates_quants$trait)]
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
  geom_abline() +
  ylab('Predicted') +
  xlab('Observed')
dev.off()


# trait fits --------------------------------------------------
pullSamples <- function(mod, param_name){
  list_of_draws <- rstan::extract(mod)
  samps <- data.frame(list_of_draws[param_name])
  sampQuantiles <- colQuantiles(list_of_draws[[param_name]], na.rm = T, probs = c(0.025, 0.50, 0.975))
  sampQuantiles <- ifelse(sampQuantiles < 0, 0, sampQuantiles)
  return(sampQuantiles)
}

plot_temperature_samples <- function(mod, param_name, df, plotPoints){
  # pull model samples
  if(grepl('R0', param_name) == FALSE){
    param_name <- paste0(param_name, '_climate_new')
  }
  x <- pullSamples(mod = mod, param_name = param_name)
  # set plotting conditions
  xvals <- df[['climate_temp_new']]
  xLabel <- expression(paste("Temperature (",degree,"C)"))
  yMax = max(x[,3])
  # plot
  plot(xvals, x[,2], type = 'l', lwd = 2, ylab = param_name, xlab = xLabel, ylim = c(0, yMax), main = param_name)
  lines(xvals,x[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(xvals, x[,3], lty=2, col='red', ylim=c(0,yMax))
  # add points
  if(plotPoints == TRUE){
    param_name_orig <- gsub('_new', '', param_name)
    param_vals <- gsub('_new', '_temp', param_name)
    points(df[[param_vals]], df[[param_name_orig]], pch=16, ylim=c(0,yMax))
  }
}


# plot temperature dependent traits and R0 model
pdf('figures/R0_Zika_trait_fit_plots.pdf', width = 11, height = 8.5)
par(mfrow = c(2, 3)) 
plot_temperature_samples(mod = r0_mod, param_name = 'alpha', df = mod_data, plotPoints = T)
plot_temperature_samples(mod = r0_mod, param_name = 'b', df = mod_data, plotPoints = T)
plot_temperature_samples(mod = r0_mod, param_name = 'pMI', df = mod_data, plotPoints = T)
plot_temperature_samples(mod = r0_mod, param_name = 'EIR', df = mod_data, plotPoints = T)
plot_temperature_samples(mod = r0_mod, param_name = 'lf', df = mod_data, plotPoints = T)
plot_temperature_samples(mod = r0_mod, param_name = 'R0_climate', df = mod_data, plotPoints = F)
dev.off()


# omega vs ancestry
aa_key <- data.frame('key' = as.character(seq(1, 66, 1)), 'aa' = mod_data$ancestry_aa_new)

list_of_draws <- rstan::extract(r0_mod)
aa_samps <- data.frame(list_of_draws['omega_ancestry_new'])
aa_samps_long <- aa_samps %>% gather()
aa_samps_long$key <- gsub('omega_ancestry_new.', '', aa_samps_long$key)
aa_samps_long <- aa_samps_long %>% left_join(aa_key) %>%
  group_by(aa) %>%
  summarise('lower' = quantile(value, 0.025)
            , 'median' = quantile(value, 0.50)
            , 'upper' = quantile(value, 0.975)
            )

pdf('figures/Omega_fit.pdf', width = 11, height = 8.5)
plot(aa_samps_long$aa, aa_samps_long$median, type = 'l', ylim= c(0,1), xlab = 'Proportion Ae. aegypti ancestry', ylab = 'Omega (prob biting humans | ancestry)')
lines(aa_samps_long$aa, aa_samps_long$lower, lty=2, col='red')
lines(aa_samps_long$aa, aa_samps_long$upper, lty=2, col='red')
points(mod_data$omega_ancestry_aa, mod_data$omega_ancestry, pch = 16)
dev.off()

# pMI ancestry modeled across strains and doses
pred_indexes <- names(r0_mod)[grep('pMI_ancestry_new', names(r0_mod))]
pred_estimates <- rstan::extract(r0_mod, pred_indexes)
pred_estimates <- lapply(pred_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)
pred_estimates_quants <- do.call(rbind.data.frame, pred_estimates)
colnames(pred_estimates_quants) <- c('lower', 'median', 'upper')

pred_estimates_quants$anc <- mod_data$ancestry_aa_new
pred_estimates_quants$Dose <- mod_data$ancestry_dose_new
pred_estimates_quants$Virus <- mod_data$ancestry_Virus_new
pred_estimates_quants$Virus <- ifelse(pred_estimates_quants$Virus == 1, 'ZIKV_Senegal_2011', 'ZIKV_Cambodia_2010')

ggplot(pred_estimates_quants, aes(x = anc, y = median, color = as.factor(Dose), fill = as.factor(Dose))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linetype = 'dashed') +
  facet_wrap(~Virus) +
  theme_classic() +
  xlab('aa ancestry') +
  ylab('Pr(Infection)')

ggsave('figures/pMI_ancestry_model.pdf')  


plot_ancestry_R0 <- function(param_name){
  pred_indexes <- names(r0_mod)[grep(param_name, names(r0_mod))]
  pred_estimates <- rstan::extract(r0_mod, pred_indexes)
  pred_estimates <- lapply(pred_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)
  pred_estimates_quants <- do.call(rbind.data.frame, pred_estimates)
  colnames(pred_estimates_quants) <- c('lower', 'median', 'upper')
  
  pred_estimates_quants$anc <- mod_data$ancestry_aa_new
  pred_estimates_quants$Dose <- mod_data$ancestry_dose_new
  pred_estimates_quants$Virus <- mod_data$ancestry_Virus_new
  pred_estimates_quants$Virus <- ifelse(pred_estimates_quants$Virus == 1, 'ZIKV_Senegal_2011', 'ZIKV_Cambodia_2010')
    
  if(param_name == 'pMI_ancestry_new'){
    yLabel = 'Pr(Infection)'
  } else {
    yLabel = 'R0'
  }
  
  ggplot(pred_estimates_quants, aes(x = anc, y = median, color = as.factor(Dose), fill = as.factor(Dose))) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linetype = 'dashed') +
    facet_wrap(~Virus) +
    theme_classic() +
    xlab('Proportion Ae. aegypti ancestry') +
    ylab(yLabel)
  
}


plot_ancestry_R0(param_name = 'pMI_ancestry_new')
ggsave('figures/pMI_ancestry_model.pdf') 

plot_ancestry_R0(param_name = 'R0_ancestry_pMI')
ggsave('figures/R0_ancestry_pMI.pdf') 

plot_ancestry_R0(param_name = 'R0_ancestry_omega')
ggsave('figures/R0_ancestry_omega.pdf') 

# need to rerun, R0 1-9 missing
plot_ancestry_R0(param_name = 'R0_ancestry') # R0_ancestry\\[.\\d]
ggsave('figures/R0_ancestry.pdf') 

# survey site R0 plots
format_survey_samples <- function(paramName){
  r0vals <- pullSamples(mod = r0_mod, param_name = paramName)
  r0vals <- as.data.frame(r0vals)
  colnames(r0vals) <- c('lower', 'median', 'upper')
  colnames(r0vals) <- paste(paramName, colnames(r0vals), sep = '_')
  r0vals$anc <- mod_data$surveys_aa
  r0vals$Virus <- mod_data$surveys_Virus
  r0vals$Dose <- mod_data$surveys_dose
  # r0vals$model <- paramName
  return(r0vals)
}

ss_r0_climate <- format_survey_samples(paramName = 'R0_climate_surveys')
ss_r0_ancestry <- format_survey_samples(paramName = 'R0_ancestry_surveys')
ss_r0_full <- format_survey_samples(paramName = 'R0_full_surveys')

ss <- ss_r0_climate %>%
  left_join(ss_r0_ancestry) %>%
  left_join(ss_r0_full)

ss$Virus <- ifelse(ss$Virus == 1, 'ZIKV_Senegal_2011', 'ZIKV_Cambodia_2010')

ggplot(ss, aes(x = R0_climate_surveys_median, y = R0_ancestry_surveys_median)) +
  geom_point() +
  # geom_errorbarh(aes(xmax = R0_climate_surveys_upper, xmin = R0_climate_surveys_lower, height = .2)) +
  # geom_errorbar(aes(ymax = R0_ancestry_surveys_upper, ymin = R0_ancestry_surveys_lower)) +
  facet_wrap(~Virus + Dose) +
  theme_bw() +
  ylim(0,4) +
  xlim(0,4) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1)

ggplot(ss, aes(x = R0_full_surveys_median, y = R0_climate_surveys_median)) +
  geom_point() +
  facet_wrap(~Virus + Dose) +
  theme_bw() +
  ylim(0,4) +
  xlim(0,4) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1)


ggplot(ss, aes(x = R0_full_surveys_median, y = R0_ancestry_surveys_median)) +
  geom_point() +
  facet_wrap(~Virus + Dose) +
  theme_bw() +
  ylim(0,4) +
  xlim(0,4) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1)


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
# pMI_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_constant', priorMean = 4.91E-04, priorSD = 0.01)
# pMI_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Tmin', priorMean = 12.22, priorSD = 20)
# pMI_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Tmax', priorMean = 37.46, priorSD = 20)
pMI_clim_rmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_rmax', priorMean = 0.24, priorSD = 0.03)
pMI_clim_Topt <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Topt', priorMean = 30.08, priorSD = 0.38)
pMI_clim_a <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_a', priorMean = 3.60, priorSD = 0.41)
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
          # , pMI_clim_const
          # , pMI_clim_Tmin
          # , pMI_clim_Tmax
          , pMI_clim_rmax
          , pMI_clim_Topt
          , pMI_clim_a
          , EIR_clim_const
          , EIR_clim_Tmin
          , EIR_clim_Tmax
          , lf_clim_const
          , lf_clim_Tmin
          , lf_clim_Tmax
          )
ggsave('figures/prior_vs_posterior_plots.pdf', p, width = 11, height = 11)


# plot survey points R0 values




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
