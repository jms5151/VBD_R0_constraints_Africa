# load libraries
library(rstan)
library(boot)
library(matrixStats)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(egg)
library(ggrepel)
library(cowplot)
library(rnaturalearth)
library(latticeExtra)
library(viridisLite)

# open model
r0_mod <- readRDS('../models/stan_model_fit_zikv.rds')
load('../VBD-data/model_data_zikv.RData')
mod_data <- model_data_zikv

# extract samples
list_of_draws <- rstan::extract(r0_mod)

# functions --------------------------------------------------------------------
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
    xLabel <- 'Proportion non-African ancestry'
    pointDataXColName <- gsub('_new', '_aa', genQuantName)
  }
  pointDataYColName <- gsub('_new', '', genQuantName)
  # set plotting conditions
  yMax = max(x[,3])
  # plot
  titleName <- gsub('_.*', '', genQuantName)
  plot(xvals, x[,2], type = 'l', lwd = 2, ylab = '', xlab = xLabel, ylim = c(0, yMax), main = titleName, las = 1)
  lines(xvals, x[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(xvals, x[,3], lty=2, col='red', ylim=c(0,yMax))
  # add points
  if(points == TRUE){
    points(mod_data[[pointDataXColName]], mod_data[[pointDataYColName]], pch=16, ylim=c(0,yMax))
  }
}

concatAncestrySamples <- function(validationName, genQuantName, percentiles){
  # get data
  x <- pullSamples(validationName = validationName, genQuantName = genQuantName, percentiles)
  x <- as.data.frame(x)
  colnames(x) <- c('lower', 'median', 'upper')
  indexes <- which(mod_data$validationtype == validationName)
  x$anc <- mod_data$aa_new[indexes]
  x$site = mod_data$location[indexes]
  x$year = mod_data$year[indexes]
  x$temp = mod_data$temp_new[indexes]
  x$lat = mod_data$lat[indexes]
  x$lon = mod_data$lon[indexes]
  return(x)  
}

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
    xlab(expression(paste('Full model ', R[0]))) +
    ylab('') +
    geom_text_repel(aes(label = site)) +
    theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm"))
}


plotSurveySites <- function(x){
  # plot
  clim_vs_full <- plotWithUncertainty(df = x, xval = 'Full_median', yval = 'Climate_median') + ylab(expression(paste('Climate model  ', R[0])))
  full_vs_anc <- plotWithUncertainty(df = x, xval = 'Full_median', yval = 'Ancestry_median') + ylab(expression(paste('Ancestry model  ', R[0])))
  # return grid of plots
  plot_grid(clim_vs_full, full_vs_anc, nrow = 1)
}

overlay_distributions_plot <- function(mod, param_name, type, priorValue1, priorValue2){
  # get posterior distribution
  df <- as.data.frame(rstan::extract(mod, param_name))
  # plot
  paramTitle <- gsub('_climate_|_ancestry_', ', ', param_name)
  p <- ggplot(df, aes_string(param_name)) +
    geom_histogram(aes(y = after_stat(density)), color = 'black', fill =  'lightblue', alpha = 0.4) +
    theme_classic() +
    ggtitle(paramTitle) +
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
            , 'pMI_ancestry_d'
            , 'pMI_ancestry_e'
            , 'pMI_ancestry_sigma'
            , 'EIR_climate_Tmin'
            , 'EIR_climate_Tmax'
            , 'EIR_climate_constant'
            , 'EIR_climate_sigma'
            , 'lf_climate_Tmin'
            , 'lf_climate_Tmax'
            , 'lf_climate_constant' 
            , 'lf_climate_sigma'
)

# create clearer labels
my_labels <- gsub('_ancestry_|_climate_', ', ', params)

# plot
plot_list <- list() 
df <- split(iris,iris$Species)

for(i in seq_along(params)){
  plot_list[[i]] <- rstan::traceplot(r0_mod, par = params[i]) + ggtitle(my_labels[i]) + ylab('') + theme(legend.position="none") + theme(plot.title = element_text(size = 12, hjust = 0.5))
}

tracePlots <- grid.arrange(grobs=plot_list, ncol = 4)

ggsave('figures/R0_stan_zikv_traceplots.pdf', tracePlots, width = 8.5, height = 11)

# ppc plots --------------------------------------------------------

# extract ppc samples and summarise
ppc_indexes <- names(r0_mod)[grep('ppc', names(r0_mod))]
ppc_estimates <- rstan::extract(r0_mod, ppc_indexes)
ppc_estimates <- lapply(ppc_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)

ppc_estimates_quants <- do.call(rbind.data.frame, ppc_estimates)
colnames(ppc_estimates_quants) <- c('lower', 'median', 'upper')
ppc_estimates_quants$trait <- names(ppc_estimates)
ppc_estimates_quants$trait <- gsub('_ppc.*', '', ppc_estimates_quants$trait)

prior_data_df <- mod_data[unique(ppc_estimates_quants$trait)]
prior_data_df <- map_df(prior_data_df, ~as.data.frame(.x), .id='trait2')
colnames(prior_data_df)[2] <- 'value'

# combine estimated and observed data
ppc <- cbind(ppc_estimates_quants, prior_data_df)
ppc$trait <- gsub('_ancestry|_climate', '', ppc$trait)
ppc$trait <- gsub('_', ', ', ppc$trait)


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
pdf('figures/R0_Zika_trait_fit_plots.pdf', width = 11, height = 7.5)
par(mfrow = c(2, 4)) 
plotParameterSamples(validationName = 'temperature', genQuantName = 'alpha_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'b_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'pMI_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'EIR_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'lf_climate_new', points = T)
plotParameterSamples(validationName = 'ancestry', genQuantName = 'omega_ancestry_new', points = T)
plotParameterSamples(validationName = 'ancestry', genQuantName = 'pMI_ancestry_new', points = T)
dev.off()

# R0 models --------------------------------------------------------------------
pdf('figures/R0_climate_ancestry.pdf', width = 8, height = 5)
par(mfrow = c(1, 2)) 
plotParameterSamples(validationName = 'temperature', genQuantName = 'R0_climate_new', points = F)
plotParameterSamples(validationName = 'ancestry', genQuantName = 'R0_ancestry_new', points = F)
dev.off()

# survey site scatterplots -----------------------------------------------------
surveys_r0_ancestry <- concatAncestrySamples(validationName = 'surveys', genQuantName = 'R0_ancestry_new', percentiles = 50)
colnames(surveys_r0_ancestry)[1:3] <- paste0('Ancestry_', colnames(surveys_r0_ancestry)[1:3])

surveys_r0_climate <- concatAncestrySamples(validationName = 'surveys', genQuantName = 'R0_climate_new', percentiles = 50)
colnames(surveys_r0_climate)[1:3] <- paste0('Climate_', colnames(surveys_r0_climate)[1:3])

surveys_r0_full <- concatAncestrySamples(validationName = 'surveys', genQuantName = 'R0_full_new', percentiles = 50)
colnames(surveys_r0_full)[1:3] <- paste0('Full_', colnames(surveys_r0_full)[1:3])

siteR0Estimates <- surveys_r0_ancestry %>% left_join(surveys_r0_climate) %>% left_join(surveys_r0_full)

sitescatter <- plotSurveySites(x = siteR0Estimates)
ggsave('figures/R0_scatterplot.pdf', sitescatter, width = 9, height = 6)

# N with R0 > 1 / model
sum(siteR0Estimates$Climate_median>1)
sum(siteR0Estimates$Ancestry_median>1)
sum(siteR0Estimates$Full_median>1)
siteR0Estimates$site[which(siteR0Estimates$Full_median>1)]

# contour plot -----------------------------------------------------------------
contour_samps <- concatAncestrySamples(validationName = 'contour', genQuantName = 'R0_full_new', percentiles = 50)

contourPlot <- ggplot(contour_samps, aes(temp, anc, z=median)) +
  geom_contour_filled() +
  guides(fill=guide_legend(expression(R[0]))) +
  metR::geom_text_contour(aes(z = median), col = 'white', size = 5) +
  theme(panel.grid=element_blank(), text=element_text(size=15)) +  # delete grid lines
  scale_x_continuous(limits=c(min(contour_samps$temp),max(contour_samps$temp)), expand=c(0,0)) + # set x limits
  scale_y_continuous(limits=c(min(contour_samps$anc),max(contour_samps$anc)), expand=c(0,0)) +  # set y limits
  xlab(expression(paste("Temperature (",degree,"C)"))) +
  ylab('Proportion non-African ancestry') +
  geom_point(surveys_r0_full, mapping = aes(x = temp, y = anc, z = 0), fill = 'black', color = 'white', pch = 16, size = 3) 

ggsave('figures/R0_countour_plot.pdf')
  
surveyValidationPlots <- ggarrange(sitescatter, contourPlot, ncol = 1)
ggsave('figures/fig1.pdf', surveyValidationPlots)

# big cities through time ------------------------------------
bc <- concatAncestrySamples(validationName = 'big_cities', genQuantName = 'R0_full_new', percentiles = 95)

bc2 <- bc[,c('site', 'year', 'median', 'lat', 'lon')] %>%
  filter(year == '1970' | year == '2020' | year == '2090-2100') %>%
  spread(key = year, value = median) 

colnames(bc2)[4:6] <- paste0('Year_', colnames(bc2)[4:6])
colnames(bc2)[6] <- 'Year_2090_2100'

# add ordered factor
bc2$Suitability <- ifelse(bc2$Year_2090_2100 > 1, 1, 0)
bc2$Suitability <- ifelse(bc2$Year_2020 > 1, 2, bc2$Suitability)
bc2$Suitability <- ifelse(bc2$Year_1970 > 1, 3, bc2$Suitability)

bc2 <- bc2[order(bc2$Suitability, bc2$Year_2090_2100),]
bc2$rank <- seq(1, nrow(bc2), 1)

bc2$site <- fct_reorder(bc2$site, bc2$rank)
bc2$Suitability <- as.factor(bc2$Suitability)

# Map
source('../google_api_key.R')

bc2$R0_max = pmax(bc2$Year_1970, bc2$Year_2020, bc2$Year_2090_2100)

world <- ne_countries(scale='medium', returnclass = 'sf')

africa <- world %>% 
  filter(continent == "Africa")

africaMap <- ggplot(data = africa) +
  geom_sf(fill = 'grey95') +
  coord_sf(xlim = c(-20, 55), ylim = c(-40, 40)) +
  geom_point(data = bc2, mapping = aes(x = lon, y = lat, color = Suitability, size = R0_max)) +
  xlab('') +
  ylab('') +
  scale_color_manual(
    values = c('navyblue', 'darkgreen', 'orange', 'maroon'),
    labels = c('Not imminent', 'Future', 'Present', 'Past'),
    guide = guide_legend(reverse = TRUE, override.aes=list(lwd = 1.3))
  ) + 
  # theme(legend.position = c(.15, .3), legend.background = element_rect(fill='transparent')) +
  theme(legend.position = 'bottom', legend.direction = 'vertical', legend.key = element_rect(fill = "transparent")) +
  scale_size_continuous(name = expression(paste('Maximum ', R[0])), breaks = seq(0,3,0.5)) +
  ### May want to change title regarding city pops
  labs(title = expression(paste('Maximum ', R[0],' between 1970 & 2100')),
       subtitle = 'Color indicates whether conditions were suitable for Zika transmission \nin the past (1970), present (2020), or future (2090-2100)') +
  scale_y_discrete(position = "right")

ggsave('figures/Africa_R0_map.pdf', africaMap, width = 7, height = 7)

# lollipop plot
bc3 <- bc %>%
  left_join(bc2[,c('site', 'Suitability', 'rank')])

bc3$site <- gsub('Democratic Republic of the Congo', 'DRC', bc3$site)
bc3$site <- fct_reorder(bc3$site, bc3$rank)

lollipopPlot <-ggplot(bc3, aes(x = median, y = site, pch = as.factor(year), group = site, color = Suitability)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  xlab(expression('R'[0])) +
  ylab('') +
  geom_vline(xintercept = 1, linetype = 2) +
  ggtitle(expression(paste('Change in ', R[0],' through time'))) +
  theme(legend.position = c(.7, .2), legend.background = element_rect(fill='transparent')) + 
  scale_color_manual(
    values = c('navyblue', 'darkgreen', 'orange', 'maroon'),
    labels = c('Not imminent', 'Future', 'Present', 'Past'),
    guide = guide_legend(reverse = TRUE, override.aes=list(linetype = 1, shape = NA, lwd = 1.3))
  ) +
  scale_shape_manual(name = 'Year', 
                     values = c('1970' = 15, '2020' = 2, '2090-2100' = 16)
  )

ggsave('figures/big_cities_lollipop_1970-2100.pdf', lollipopPlot, width = 11, height = 11)

# combine
bigCities <- plot_grid(lollipopPlot, africaMap, ncol = 2)

# save
ggsave('figures/big_cities_1970-2100.pdf', bigCities, height = 7.5, width = 11)

# prior vs posterior plots -----------------------------------
a_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_constant', type = 'normal', priorValue1 = 2.02E-04, priorValue2 = 0.01)
a_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmin', type = 'normal', priorValue1 =  13.35, priorValue2 = 20)
a_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmax', type = 'normal', priorValue1 = 40.08, priorValue2 = 20)
a_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
b_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_constant', type = 'normal', priorValue1 = 8.49E-04, priorValue2 = 0.01)
b_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmin', type = 'normal', priorValue1 = 17.05, priorValue2 = 20)
b_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmax', type = 'normal', priorValue1 = 35.83, priorValue2 = 20)
b_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
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
