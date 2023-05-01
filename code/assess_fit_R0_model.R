# load libraries
library(rstan)
# library(shinystan)
library(boot)
library(matrixStats)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(cowplot)

# open model
r0_mod <- readRDS('../models/stan_model_fit_zikv.rds')
load('../VBD-data/model_data_zikv.RData')
mod_data <- model_data_zikv

# assess model fit
# launch_shinystan(r0_mod)

# trace plots ------------------------------------------------------
# list parameters
# ignores deltas, etc.
params <- c('omega_ancestry_constant'
            , 'omega_ancestry_d'
            , 'omega_ancestry_e'
            , 'omega_ancestry_sigma'
            , 'alpha_climate_Tmin'
            , 'alpha_climate_Tmax'
            , 'alpha_climate_constant'
            , 'alpha_climate_sigma'
            , 'b_climate_Tmin'
            , 'b_climate_Tmax'
            , 'b_climate_constant'
            , 'b_climate_sigma'
            # , 'pMI_climate_Tmin'
            # , 'pMI_climate_Tmax'
            # , 'pMI_climate_constant'
            # , 'pMI_climate_sigma'
            , 'pMI_climate_rmax'
            , 'pMI_climate_Topt'
            , 'pMI_climate_a'
            , 'pMI_climate_sigma'
            , 'pMI_ancestry_b0'
            , 'pMI_ancestry_beta[1]'
            , 'pMI_ancestry_beta[2]'
            , 'pMI_ancestry_gamma'
            , 'EIR_climate_Tmin'
            , 'EIR_climate_Tmax'
            , 'EIR_climate_constant'
            , 'EIR_climate_sigma'
            , 'lf_climate_Tmin'
            , 'lf_climate_Tmax'
            , 'lf_climate_constant' 
            , 'lf_climate_sigma'
)
pdf('figures/R0_stan_zikv_traceplots.pdf', width = 11, height = 8.5)
rstan::traceplot(r0_mod, par = c('lp__', params), ncol = 5, nrow = 6)
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
list_of_draws <- rstan::extract(r0_mod)

pullSamples <- function(validationName, genQuantName, percentiles = 95){
  indexes <- which(mod_data$validationtype == validationName)
  samps <- list_of_draws[[genQuantName]][, indexes]
  if(percentiles == 50){
    sampQuantiles <- colQuantiles(samps, na.rm = T, probs = c(0.25, 0.50, 0.75))
  } else {
    sampQuantiles <- colQuantiles(samps, na.rm = T, probs = c(0.025, 0.50, 0.975))
  }
  sampQuantiles <- ifelse(sampQuantiles < 0, 0, sampQuantiles)
  return(sampQuantiles)
}

plotParameterSamples <- function(validationName, genQuantName, points){
  # get samples
  x <- pullSamples(validationName = validationName, genQuantName = genQuantName)
  # get data
  indexes <- which(mod_data$validationtype == validationName)
  if(grepl('climate', genQuantName) == TRUE){
    xvals <- mod_data$temp_new[indexes]
    xLabel <- expression(paste("Temperature (",degree,"C)"))
    pointDataXColName <- gsub('_new', '_temp', genQuantName)
  } else if(grepl('ancestry', genQuantName) == TRUE){
    xvals <- mod_data$aa_new[indexes]
    xLabel <- 'Proportion Ae. aegypti ancestry'
    pointDataXColName <- gsub('_new', '_aa', genQuantName)
  }
  pointDataYColName <- gsub('_new', '', genQuantName)
  # set plotting conditions
  yMax = max(x[,3])
  # plot
  plot(xvals, x[,2], type = 'l', lwd = 2, ylab = genQuantName, xlab = xLabel, ylim = c(0, yMax), main = genQuantName)
  lines(xvals, x[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(xvals, x[,3], lty=2, col='red', ylim=c(0,yMax))
  # add points
  if(points == TRUE){
    points(mod_data[[pointDataXColName]], mod_data[[pointDataYColName]], pch=16, ylim=c(0,yMax))
  }
}

pdf('figures/R0_Zika_trait_fit_plots.pdf', width = 11, height = 8.5)
par(mfrow = c(2, 3)) 
plotParameterSamples(validationName = 'temperature', genQuantName = 'alpha_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'b_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'pMI_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'EIR_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'lf_climate_new', points = T)
plotParameterSamples(validationName = 'ancestry', genQuantName = 'omega_ancestry_new', points = T)
dev.off()

# pMI ancestry modeled across strains and doses ---------------------------------

concatAncestrySamples <- function(validationName, genQuantName, percentiles){
  # get data
  x <- pullSamples(validationName = validationName, genQuantName = genQuantName, percentiles)
  x <- as.data.frame(x)
  colnames(x) <- c('lower', 'median', 'upper')
  indexes <- which(mod_data$validationtype == validationName)
  x$anc <- mod_data$aa_new[indexes]
  x$Dose <- mod_data$dose_new[indexes]
  x$Virus <- mod_data$Virus_new[indexes]
  x$Virus <- ifelse(x$Virus == 1, 'ZIKV_Senegal_2011', 'ZIKV_Cambodia_2010')
  return(x)  
}

ancestryFitsAcrossTreatments <- function(validationName, genQuantName){
  # get data
  x <- concatAncestrySamples(validationName = validationName, genQuantName = genQuantName, percentiles = 95)
  # plot labels
  if(genQuantName == 'pMI_ancestry_new'){
    yLabel = 'Pr(Infection)'
  } else {
    yLabel = 'R0'
  }
  # plot
  ggplot(x, aes(x = anc, y = median, color = as.factor(Dose), fill = as.factor(Dose))) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linetype = 'dashed') +
    facet_wrap(~Virus) +
    theme_classic() +
    xlab('aa ancestry') +
    ylab(yLabel)
  
}

ancestryFitsAcrossTreatments(validationName = 'ancestry', genQuantName = 'pMI_ancestry_new')
ggsave('figures/pMI_ancestry_model.pdf')  


# R0 models --------------------------------------------------------------------
ancestryFitsAcrossTreatments(validationName = 'ancestry', genQuantName = 'R0_ancestry_new')
ggsave('figures/R0_ancestry.pdf')

pdf('figures/R0_climate.pdf', width = 8, height = 6)
plotParameterSamples(validationName = 'temperature', genQuantName = 'R0_climate_new', points = F)
dev.off()

# survey site scatterplots -----------------------------------------------------
surveys_r0_ancestry <- concatAncestrySamples(validationName = 'surveys', genQuantName = 'R0_ancestry_new', percentiles = 50)
# surveys_r0_ancestry$model <- 'ancestry'
colnames(surveys_r0_ancestry) <- paste0('Ancestry_', colnames(surveys_r0_ancestry))

surveys_r0_climate <- concatAncestrySamples(validationName = 'surveys', genQuantName = 'R0_climate_new', percentiles = 50)
# surveys_r0_climate$model <- 'climate'
colnames(surveys_r0_climate) <- paste0('Climate_', colnames(surveys_r0_climate))

surveys_r0_full <- concatAncestrySamples(validationName = 'surveys', genQuantName = 'R0_full_new', percentiles = 50)
# surveys_r0_full$model <- 'full'
colnames(surveys_r0_full) <- paste0('Full_', colnames(surveys_r0_full))

# siteR0Estimates <- do.call(rbind, list(surveys_r0_ancestry, surveys_r0_climate, surveys_r0_full))
siteR0Estimates <- do.call(cbind, list(surveys_r0_ancestry, surveys_r0_climate, surveys_r0_full))


############### STOPPED HERE
############### NEED TO FIGURE OUT HOW TO DO 3-PLOT
gridPlotDoseVirus <- function(df, xval, yval){
  ggplot(df, aes_string(x = df[,xval], y = df[,yval])) +
    geom_point() +
    facet_wrap(~Virus + Dose) +
    theme_bw() +
    ylim(0,4) +
    xlim(0,4) +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = 1) +
    xlab(xval) +
    ylab(yval) 
}


plotSurveySites <- function(virus, dose){
    
}


df <- subset(siteR0Estimates, Virus == 'ZIKV_Cambodia_2010' & Dose == 1275000)

ggplot(df, aes())
# siteR0Estimates <- do.call(cbind, list(surveys_r0_ancestry, surveys_r0_climate, surveys_r0_full))

# format_survey_samples <- function(paramName){
#   r0vals <- pullSamples50(mod = r0_mod, param_name = paramName)
#   r0vals <- as.data.frame(r0vals)
#   colnames(r0vals) <- c('lower', 'median', 'upper')
#   colnames(r0vals) <- paste(paramName, colnames(r0vals), sep = '_')
#   r0vals$anc <- mod_data$surveys_aa
#   r0vals$Virus <- mod_data$surveys_Virus
#   r0vals$Dose <- mod_data$surveys_dose
#   r0vals$site <- mod_data$surveys_location
#   # r0vals$model <- paramName
#   return(r0vals)
# }
# 
# ss_r0_climate <- format_survey_samples(paramName = 'R0_climate_surveys')
# ss_r0_ancestry <- format_survey_samples(paramName = 'R0_ancestry_surveys')
# ss_r0_full <- format_survey_samples(paramName = 'R0_full_surveys')
# 
# ss <- ss_r0_climate %>%
#   left_join(ss_r0_ancestry) %>%
#   left_join(ss_r0_full)

# ss$Virus <- ifelse(ss$Virus == 1, 'ZIKV_Senegal_2011', 'ZIKV_Cambodia_2010')

gridPlotDoseVirus <- function(df, xval, yval){
  ggplot(df, aes_string(x = df[,xval], y = df[,yval])) +
    geom_point() +
    facet_wrap(~Virus + Dose) +
    theme_bw() +
    ylim(0,4) +
    xlim(0,4) +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = 1) +
    xlab(xval) +
    ylab(yval) 
}

gridPlotDoseVirus(df = cambodia_low_dose, xval = 'R0_climate_surveys_median', yval = 'R0_ancestry_surveys_median')
ggsave('figures/R0_surveys_climate_vs_ancestry.pdf', height = 8.5, width = 11) 

gridPlotDoseVirus(df = ss, xval = 'R0_full_surveys_median', yval = 'R0_climate_surveys_median')
ggsave('figures/R0_surveys_full_vs_climate.pdf', height = 8.5, width = 11) 

gridPlotDoseVirus(df = ss, xval = 'R0_full_surveys_median', yval = 'R0_ancestry_surveys_median')
ggsave('figures/R0_surveys_full_vs_ancestry.pdf', height = 8.5, width = 11) 

plotWithUncertainty <- function(df, xval, yval){
  xupper <- gsub('_median', '_upper', xval)
  xlower <- gsub('_median', '_lower', xval)
  yupper <- gsub('_median', '_upper', yval)
  ylower <- gsub('_median', '_lower', yval)
  df$site <- ifelse(df[,xval] > 1 & df[,yval] > 1, df$site, '')

  ggplot(df, aes_string(x = df[,xval], y = df[,yval])) +
    geom_errorbarh(aes_string(xmax = df[,xupper], xmin = df[,xlower]), col = 'darkgrey') +
    geom_errorbar(aes_string(ymax = df[,yupper], ymin = df[,ylower]), col = 'darkgrey') +
    geom_point(size = 2, color = 'black') +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    ylim(0,6) +
    xlim(0,6) +
    geom_hline(yintercept = 1,linetype=2) +
    geom_vline(xintercept = 1,linetype=2) +
    xlab(gsub('_|surveys|median', ' ', xval)) +
    ylab(gsub('_|surveys|median', ' ', yval)) +
    geom_text_repel(aes(label = site))
}

# highlight since dose and site
cambodia_low_dose <- subset(ss, Virus == 'ZIKV_Cambodia_2010' & Dose == 1275000)

cam1 <- plotWithUncertainty(df = cambodia_low_dose, xval = 'R0_climate_surveys_median', yval = 'R0_ancestry_surveys_median')
# ggsave('figures/R0_surveys_climate_vs_ancestry_cambodia_low_dose.pdf', height = 8.5, width = 11) 

cam2 <- plotWithUncertainty(df = cambodia_low_dose, xval = 'R0_full_surveys_median', yval = 'R0_climate_surveys_median')
# ggsave('figures/R0_surveys_full_vs_climate_cambodia_low_dose.pdf', height = 8.5, width = 11) 

cam3 <- plotWithUncertainty(df = cambodia_low_dose, xval = 'R0_full_surveys_median', yval = 'R0_ancestry_surveys_median')
# ggsave('figures/R0_surveys_full_vs_ancestry_cambodia_low_dose.pdf', height = 8.5, width = 11) 

cam <- plot_grid(cam1, cam2, cam3, nrow = 1)
ggsave('figures/R0_scatterplot_cambodia_low_dose.pdf.pdf', cam, width = 12, height = 4)

# summary
cambodia_low_dose_summary <- cambodia_low_dose %>%
  summarise(climate_N = sum(R0_climate_surveys_median > 1)
            , ancestry_N = sum(R0_ancestry_surveys_median > 1)
            , full_N = sum(R0_full_surveys_median > 1))

# prior vs posterior plots -----------------------------------

# function to overlay histogram and normal density
overlay_distributions_plot <- function(mod, param_name, type, priorValue1, priorValue2){
  # get posterior distribution
  df <- as.data.frame(rstan::extract(mod, param_name))
  # plot
  p <- ggplot(df, aes_string(param_name)) +
    geom_histogram(aes(y = after_stat(density)), color = 'black', fill =  'yellow', alpha = 0.4) +
    theme_classic() +
    ggtitle(param_name) +
    ylab('') +
    xlab('')
  
  if(type == 'normal'){
    p  + 
      stat_function(
        fun = dnorm, 
        args = list(mean = priorValue1, sd = priorValue2), 
        lwd = 2, 
        col = 'black'
      )
  } else if(type == 'uniform'){
    p  + 
      stat_function(
        fun = dunif, 
        args = list(min = priorValue1, max = priorValue2), 
        lwd = 2, 
        col = 'black'
      )
  }
}

# overlay_distributions_plot <- function(mod, param_name, priorMean, priorSD){
#   # get posterior distribution
#   df <- as.data.frame(rstan::extract(mod, param_name))
#   # plot
#   ggplot(df, aes_string(param_name)) +
#     geom_histogram(aes(y = after_stat(density)), color = 'black', fill =  'yellow', alpha = 0.4) +
#     stat_function(
#       fun = dunif, 
#       args = list(min = priorMin, max = priorMax), 
#       lwd = 2, 
#       col = 'black'
#     ) +
#     theme_classic() +
#     ggtitle(param_name) +
#     ylab('') +
#     xlab('')
# }

a_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_constant', type = 'normal', priorValue1 = 2.02E-04, priorValue2 = 0.01)
a_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmin', type = 'normal', priorValue1 =  13.35, priorValue2 = 20)
a_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmax', type = 'normal', priorValue1 = 40.08, priorValue2 = 20)
a_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
b_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_constant', type = 'normal', priorValue1 = 8.49E-04, priorValue2 = 0.01)
b_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmin', type = 'normal', priorValue1 = 17.05, priorValue2 = 20)
b_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmax', type = 'normal', priorValue1 = 35.83, priorValue2 = 20)
b_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
# pMI_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_constant', type = 'normal', priorValue1 = 4.91E-04, priorValue2 = 0.01)
# pMI_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Tmin', type = 'normal', priorValue1 = 12.22, priorValue2 = 20)
# pMI_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Tmax', type = 'normal', priorValue1 = 37.46, priorValue2 = 20)
pMI_clim_rmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_rmax', type = 'normal', priorValue1 = 0.24, priorValue2 = 0.03)
pMI_clim_Topt <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Topt', type = 'normal', priorValue1 = 30.08, priorValue2 = 0.38)
pMI_clim_a <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_a', type = 'normal', priorValue1 = 3.60, priorValue2 = 0.41)
pMI_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
EIR_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_constant', type = 'normal', priorValue1 = 6.65E-05, priorValue2 = 0.01)
EIR_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_Tmin', type = 'normal', priorValue1 = 10.68, priorValue2 = 20)
EIR_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_Tmax', type = 'normal', priorValue1 = 45.90, priorValue2 = 20)
EIR_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
lf_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_constant', type = 'normal', priorValue1 = -1.48E-01, priorValue2 = 0.1)
lf_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_Tmin', type = 'normal', priorValue1 = 9.16, priorValue2 = 20)
lf_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_Tmax', type = 'normal', priorValue1 = 37.73, priorValue2 = 20)
lf_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)

p <- plot_grid(a_clim_const
          , a_clim_Tmin
          , a_clim_Tmax
          , a_clim_sig
          , b_clim_const
          , b_clim_Tmin
          , b_clim_Tmax
          , b_clim_sig
          # , pMI_clim_const
          # , pMI_clim_Tmin
          # , pMI_clim_Tmax
          , pMI_clim_rmax
          , pMI_clim_Topt
          , pMI_clim_a
          , pMI_clim_sig
          , EIR_clim_const
          , EIR_clim_Tmin
          , EIR_clim_Tmax
          , EIR_clim_sig
          , lf_clim_const
          , lf_clim_Tmin
          , lf_clim_Tmax
          , lf_clim_sig
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
