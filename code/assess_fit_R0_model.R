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
library(ggpubr)

# open model
r0_mod <- readRDS('../models/stan_model_fit_zikv.rds')
load('../VBD-data/model_data_zikv.RData')
mod_data <- model_data_zikv

# summary statistics ---------------------------------------------
# Extract summary from the STAN fit object
fit_summary <- summary(r0_mod)

# Extract the summary statistics
summary_stats <- fit_summary$summary
summary_stats <- summary_stats[1:28, ]

# Extract Rhat and ESS
rhat <- round(summary_stats[, "Rhat"], 2)
ess <- round(summary_stats[, "n_eff"])

# Create a dataframe with the diagnostics
diagnostics_df <- data.frame(Parameter = rownames(summary_stats), Rhat = rhat, ESS = ess)

# save
write.csv(diagnostics_df, '../VBD-data/model_diagnostics.csv', row.names = F)

# extract samples ------------------------------------------------
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

r0_plot <- function(validationName, genQuantName) {
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
    xLabel <- 'Prop. non-African ancestry'
    pointDataXColName <- gsub('_new', '_aa', genQuantName)
  }
  p <- ggplot() +
    geom_line(aes(xvals, x[,2])) +
    geom_line(aes(xvals, x[,1]), col = 'red', linetype = 'dashed') +
    geom_line(aes(xvals, x[,3]), col = 'red', linetype = 'dashed') +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    ylab(expression(R[0])) +
    xlab(xLabel)
  
  return(p)
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
    ylim(0,4) +
    xlim(0,4) +
    geom_hline(yintercept = 1,linetype=2) +
    geom_vline(xintercept = 1,linetype=2) +
    geom_text_repel(aes(label = site)) +
    theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm"))
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

validation_plot <- function(df, xval, yval){
  xupper <- 'sePrev_Upper'
  xlower <- 'sePrev_Lower'
  yupper <- gsub('_median', '_upper', yval)
  ylower <- gsub('_median', '_lower', yval)
  
  # Calculate regression line and R-squared value
  fit <- lm(df[, yval] ~ df[, xval])
  rsquared <- summary(fit)$r.squared
  
  ggplot(df, aes_string(x = xval, y = yval)) +
    # geom_errorbarh(aes(xmax = df[, xupper], xmin = df[, xlower]), col = 'darkgrey') +
    geom_errorbar(aes(ymax = df[, yupper], ymin = df[, ylower]), col = 'darkgrey') +
    geom_point(size = 2, color = 'black') +
    geom_smooth(method = 'lm', se = FALSE, color = 'black', fullrange = TRUE) + # Extend the line
    annotate("text", x = 2, y = 3.4, # Move text to upper left
             label = bquote(paste("R"^2, " = ", .(sprintf("%.2f", rsquared)))), 
             hjust = 0, vjust = 1, col = 'black') +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    xlab('Seroprevalence') +
    ylab(expression(paste('Model estimated  ', R[0]))) +
    ylim(0,3.5) +
    geom_text_repel(aes(label = site)) +
    theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm")) +
    ggtitle(gsub('_median', ' model', yval))
}


# trace plots ------------------------------------------------------
# list parameters
# ignores deltas, etc.
params <- c('omega_ancestry_constant'
            , 'omega_ancestry_d'
            , 'omega_ancestry_e'
            # , 'omega_ancestry_sigma'
            , 'alpha_climate_Tmin'
            , 'alpha_climate_Tmax'
            , 'alpha_climate_constant'
            # , 'alpha_climate_sigma'
            , 'b_climate_Tmin'
            , 'b_climate_Tmax'
            , 'b_climate_constant'
            # , 'b_climate_sigma'
            , 'EIR_climate_Tmin'
            , 'EIR_climate_Tmax'
            , 'EIR_climate_constant'
            # , 'EIR_climate_sigma'
            , 'lf_climate_Tmin'
            , 'lf_climate_Tmax'
            , 'lf_climate_constant' 
            # , 'lf_climate_sigma'
            , 'pMI_ancestry_constant'
            , 'pMI_ancestry_d'
            , 'pMI_ancestry_e'
            # , 'pMI_ancestry_sigma'
            , 'pMI_climate_rmax'
            , 'pMI_climate_Topt'
            , 'pMI_climate_a'
)

# create clearer labels
my_labels <- gsub('_ancestry_|_climate_', ', ', params)

# plot
plot_list <- list() 

for(i in seq_along(params)){
  plot_list[[i]] <- rstan::traceplot(r0_mod, par = params[i]) + ggtitle(my_labels[i]) + ylab('') + theme(legend.position="none") + theme(plot.title = element_text(size = 12, hjust = 0.5))
}

tracePlots <- grid.arrange(grobs=plot_list, ncol = 4)

ggsave('figures/R0_stan_zikv_traceplots.pdf', tracePlots, width = 8.5, height = 11)

# list 95% CI for each parameter
param95ci <- data.frame('param_name' = character(), 'ci5' = numeric(), 'ci50' = numeric(), 'ci95' = numeric())

for(i in 1:length(params)){
  cis <- quantile(list_of_draws[[params[i]]], c(0.05, 0.50, 0.95))
  param95ci <- param95ci %>% 
    add_row(param_name = params[i], 
            'ci5' = unname(cis[1]), 
            'ci50' = unname(cis[2]),
            'ci95' = unname(cis[3])
    )
}

write.csv(param95ci, '../trait_fits.csv', row.names = F)

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

# edit pMI since there are two 
ppc$trait[ppc$trait2 == 'pMI_climate'] <- 'pMI (climate)'
ppc$trait[ppc$trait2 == 'pMI_ancestry'] <- 'pMI (ancestry)'

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
plotParameterSamples(validationName = 'ancestry', genQuantName = 'omega_ancestry_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'alpha_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'b_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'EIR_climate_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'lf_climate_new', points = T)
plotParameterSamples(validationName = 'ancestry', genQuantName = 'pMI_ancestry_new', points = T)
plotParameterSamples(validationName = 'temperature', genQuantName = 'pMI_climate_new', points = T)
dev.off()

# R0 models --------------------------------------------------------------------
# pdf('figures/R0_climate_ancestry.pdf', width = 8, height = 5)
# par(mfrow = c(1, 2)) 
# plotParameterSamples(validationName = 'ancestry', genQuantName = 'R0_ancestry_new', points = F)
# plotParameterSamples(validationName = 'temperature', genQuantName = 'R0_climate_new', points = F)
# dev.off()
anc_r0_plot <- r0_plot(validationName = 'ancestry', genQuantName = 'R0_ancestry_new') + ylim(0,6)
temp_r0_plot <- r0_plot(validationName = 'temperature', genQuantName = 'R0_climate_new') + ylim(0,6)

# MPP <- read.csv('../VBD-data/seroSites_Aaa.csv')
# MPP$Match_Site <- paste(MPP$Site, MPP$Country, sep = ', ')
# 
# MPP$Ancestry <- ifelse(MPP$aaa2015 <= 0.001, 'Low', NA)
# MPP$Ancestry <- ifelse(MPP$aaa2015 > 0.001 & MPP$aaa2015 <= 0.01, 'Transitioning', MPP$Ancestry)
# MPP$Ancestry <- ifelse(MPP$aaa2015 > 0.01 & MPP$aaa2015 <= 0.05, 'Moderate', MPP$Ancestry)
# # MPP$Ancestry <- ifelse(MPP$aaa2015 > 0.05 & MPP$aaa2015 <= 0.15, 'High', MPP$Ancestry)
# MPP$Ancestry <- ifelse(MPP$aaa2015 > 0.05, 'Strong specialization', MPP$Ancestry)
# 
# MPP$Ancestry <- factor(MPP$Ancestry, levels = c('Low', 'Transitioning', 'Moderate', 'Strong specialization'))
# MPP$Physicians_per_100000[is.na(MPP$Physicians_per_100000)] <-  min(MPP$Physicians_per_100000, na.rm = T)
# m <- MPP
# m <- subset(MPP, Neutralizing_antibodies == 'Yes')
# m <- subset(m, Country != 'Madagascar')
# m <- subset(m, Site != 'Oyo State')
# m$phys <- m$Physicians_per_100000/5
# 
# # points(v$bio8_20, v$Seroprevalence/20, pch = 16, ylim = c(0,100), xlim = c(10, 40))
# sero_plot_T <- temp_r0_plot + 
#   geom_point(aes(x = m$bio8_20, y  = m$Seroprevalence/8, color = m$Ancestry), size = m$phys) +
#   scale_y_continuous(sec.axis = sec_axis(~ .*8, name = "Seroprevalence")) +
#   labs(color = 'Non-African\nmosquito ancestry')
# 
# ggplot() +
#   geom_point(aes(x = m$bio8_20, y  = m$Seroprevalence/8, color = m$Ancestry), size = m$phys) +
#   labs(color = 'Non-African\nmosquito ancestry') +
#   theme_bw()
# 
# 
# ggsave('figures/sero_v_r0_temp.pdf', plot = sero_plot_T, width = 10, height = 6)
# 
# 
# m$Tsuit <- ifelse(m$bio8_20 <=21 | m$bio8_20 > 31, 'Low', NA)
# m$Tsuit <- ifelse(m$bio8_20 > 21 & m$bio8_20 <= 26, 'Moderate', m$Tsuit)
# m$Tsuit <- ifelse(m$bio8_20 > 26 & m$bio8_20 <= 31, 'High', m$Tsuit)
# 
# sero_plot_A <- anc_r0_plot +
#   geom_point(aes(x = m$aaa2015, y = m$Seroprevalence/8, color = m$Tsuit), size = m$phys) +
#   labs(color = 'Temperature suitability')
# 
# ggsave('figures/sero_v_r0_anc.pdf', plot = sero_plot_A, width = 10, height = 6)

# contour plot -------------
contour_samps <- concatAncestrySamples(validationName = 'contour', genQuantName = 'R0_full_new', percentiles = 50)

# With seroprevalence surveys
contourPlot <- ggplot(contour_samps, aes(temp, anc, z=median)) +
  geom_contour_filled(breaks = c(0, 0.5, 0.9, 1, 1.1, 1.5, 2, 3, 4)) + #seq(from = 0, to = 5, by = 0.5)
  guides(fill=guide_legend(expression(R[0]))) +
  metR::geom_text_contour(aes(z = median),  col = 'white', size = 5) +
  theme(panel.grid=element_blank(), text=element_text(size=15)) +  # delete grid lines
  scale_x_continuous(limits=c(min(contour_samps$temp),max(contour_samps$temp)), expand=c(0,0)) + # set x limits
  scale_y_continuous(limits=c(min(contour_samps$anc),max(contour_samps$anc)), expand=c(0,0)) +  # set y limits
  xlab(expression(paste("Temperature (",degree,"C)"))) +
  ylab('Proportion non-African ancestry') +
  geom_point(MPP, mapping = aes(x = bio8_20, y = aaa2015, z = 0), fill = 'black', color = 'white', pch = 16, size = 3) #+


r0_contours <- plot_grid(anc_r0_plot, temp_r0_plot, contourPlot, ncol = 3, rel_widths = c(0.5, 0.5, 0.9))
ggsave('figures/r0_and_contours.pdf', r0_contours, width = 12, height = 4)


# combine scatter and contour plots and save
# surveyValidationPlots <- ggarrange(Full_v_Anc, Full_v_Clim, Omega_v_pMI, contourPlot, ncol = 2, nrow = 2, labels = c('A', 'B', 'C', 'D'))
# ggsave('figures/fig1.pdf', surveyValidationPlots, width = 12, height = 9)

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
  coord_sf(xlim = c(-20, 55), ylim = c(-50, 50)) +
  geom_point(data = bc2, mapping = aes(x = lon, y = lat, color = Suitability, size = R0_max)) +
  xlab('') +
  ylab('') +
  scale_color_manual(
    values = c('navyblue', 'darkgreen', 'orange', 'maroon'),
    labels = c('Not imminent', 'Future', 'Present', 'Past'),
    guide = guide_legend(reverse = TRUE, override.aes=list(lwd = 1.3))
  ) + 
  theme_bw() +
  theme(legend.position = c(.15, .3), legend.background = element_rect(fill='transparent')) +
  # theme(legend.position = 'bottom', legend.direction = 'vertical', legend.key = element_rect(fill = "transparent")) +
  scale_size_continuous(name = expression(paste('Maximum ', R[0])), breaks = seq(0,3,0.5)) +
  ### May want to change title regarding city pops
  labs(title = expression(paste('Maximum ', R[0],' between 1970 & 2100')),
       # subtitle = 'Color indicates whether conditions were suitable for Zika transmission \nin the past (1970), present (2020), or future (2090-2100)'
       ) +
  scale_y_discrete(position = "right") 


# lollipop plot
bc3 <- bc %>%
  left_join(bc2[,c('site', 'Suitability', 'rank')])

bc3$site <- gsub('Democratic Republic of the Congo', 'DRC', bc3$site)
bc3$site <- fct_reorder(bc3$site, bc3$rank)

lollipopPlot <- ggplot(bc3, aes(x = median, y = site, pch = as.factor(year), group = site, color = Suitability)) +
  geom_point(size = 2) +
  geom_line() +
  theme_bw() +
  xlab(expression('R'[0])) +
  ylab('') +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(title = expression(paste('Change in ', R[0],' through time'))) +
  theme(legend.position = c(.8, .2), legend.background = element_rect(fill='transparent')) + 
  scale_color_manual(
    values = c('navyblue', 'darkgreen', 'orange', 'maroon'),
    labels = c('Not imminent', 'Future', 'Present', 'Past'),
    guide = guide_legend(reverse = TRUE, override.aes=list(linetype = 1, shape = NA, lwd = 1.3))
  ) +
  scale_shape_manual(name = 'Year', 
                     values = c('1970' = 15, '2020' = 2, '2090-2100' = 16)
  )

# combine
bigCities <- ggarrange(lollipopPlot, africaMap, ncol = 2, labels = c('A', 'B'))
# bigCities <- plot_grid(lollipopPlot, africaMap, ncol = 2)

# save
ggsave('figures/big_cities_1970-2100.pdf', bigCities, height = 7.5, width = 14)


ggplot(bc3, aes(x = median, y = site, pch = as.factor(year), group = site, color = Suitability)) +
  geom_point(size = 2) +
  # geom_errorbarh(aes(xmin = S_5th, xmax = S_95th), col = 'darkgrey') +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = Suitability)) +
  facet_grid(~as.factor(year)) +
  theme_bw() +
  xlab(expression('R'[0])) +
  ylab('') +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(title = expression(paste('Change in ', R[0],' through time'))) +
  theme(legend.position = c(.8, .2), legend.background = element_rect(fill='transparent')) + 
  scale_color_manual(
    values = c('navyblue', 'darkgreen', 'orange', 'maroon'),
    labels = c('Not imminent', 'Future', 'Present', 'Past'),
    guide = guide_legend(reverse = TRUE, override.aes=list(linetype = 1, shape = NA, lwd = 1.3))
  ) +
  scale_shape_manual(name = 'Year', 
                     values = c('1970' = 15, '2020' = 2, '2090-2100' = 16)
  )



# prior vs posterior plots -----------------------------------
o_anc_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'omega_ancestry_constant', type = 'normal', priorValue1 = 0.5, priorValue2 = 1)
o_anc_d <- overlay_distributions_plot(mod = r0_mod, param_name = 'omega_ancestry_d', type = 'normal', priorValue1 = 0.5, priorValue2 = 1)
o_anc_e <- overlay_distributions_plot(mod = r0_mod, param_name = 'omega_ancestry_e', type = 'normal', priorValue1 = 0.5, priorValue2 = 1)
# o_anc_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'omega_ancestry_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
a_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_constant', type = 'normal', priorValue1 = 2.02E-04, priorValue2 = 0.01)
a_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmin', type = 'normal', priorValue1 =  13.35, priorValue2 = 1)
a_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_Tmax', type = 'normal', priorValue1 = 40.08, priorValue2 = 1)
# a_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'alpha_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
b_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_constant', type = 'normal', priorValue1 = 8.49E-04, priorValue2 = 1)
b_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmin', type = 'normal', priorValue1 = 17.05, priorValue2 = 1)
b_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_Tmax', type = 'normal', priorValue1 = 35.83, priorValue2 = 1)
# b_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'b_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
EIR_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_constant', type = 'normal', priorValue1 = 6.65E-05, priorValue2 = 1)
EIR_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_Tmin', type = 'normal', priorValue1 = 10.68, priorValue2 = 1)
EIR_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_Tmax', type = 'normal', priorValue1 = 45.90, priorValue2 = 1)
# EIR_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'EIR_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
lf_clim_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_constant', type = 'normal', priorValue1 = -1.48E-01, priorValue2 = 1)
lf_clim_Tmin <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_Tmin', type = 'normal', priorValue1 = 9.16, priorValue2 = 1)
lf_clim_Tmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_Tmax', type = 'normal', priorValue1 = 37.73, priorValue2 = 1)
# lf_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'lf_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
pMI_anc_const <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_ancestry_constant', type = 'normal', priorValue1 = 0.5, priorValue2 = 1)
pMI_anc_d <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_ancestry_d', type = 'normal', priorValue1 = 0.5, priorValue2 = 1)
pMI_anc_e <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_ancestry_e', type = 'normal', priorValue1 = 0.5, priorValue2 = 1)
# pMI_anc_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_ancestry_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)
pMI_clim_rmax <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_rmax', type = 'normal', priorValue1 = 0.24, priorValue2 = 1)
pMI_clim_Topt <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_Topt', type = 'normal', priorValue1 = 30.08, priorValue2 = 1)
pMI_clim_a <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_a', type = 'normal', priorValue1 = 3.60, priorValue2 = 1)
# pMI_clim_sig <- overlay_distributions_plot(mod = r0_mod, param_name = 'pMI_climate_sigma', type = 'uniform', priorValue1 = 0, priorValue2 = 100)

p <- plot_grid(o_anc_const
               , o_anc_d
               , o_anc_e
               # , o_anc_sig
               , a_clim_const
               , a_clim_Tmin
               , a_clim_Tmax
               # , a_clim_sig
               , b_clim_const
               , b_clim_Tmin
               , b_clim_Tmax
               # , b_clim_sig
               , EIR_clim_const
               , EIR_clim_Tmin
               , EIR_clim_Tmax
               # , EIR_clim_sig
               , lf_clim_const
               , lf_clim_Tmin
               , lf_clim_Tmax
               # , lf_clim_sig
               , pMI_anc_const
               , pMI_anc_d
               , pMI_anc_e
               # , pMI_anc_sig
               , pMI_clim_rmax
               , pMI_clim_Topt
               , pMI_clim_a
               # , pMI_clim_sig
)
ggsave('figures/prior_vs_posterior_plots.pdf', p, width = 11, height = 11)

# seroprevalence validation ---------------------------------------------------
MPP <- read.csv('../VBD-data/seroSites_Aaa_v2.csv')

create_summary_dataset <- function(val_type){
  # format names
  sero_full <- concatAncestrySamples(validationName = 'seroprevalence', genQuantName = val_type, percentiles = 95)
  sero_full$Country <- MPP$Country
  sero_full$Seroprevalence <- MPP$Seroprevalence
  sero_full$NeutAnti <- MPP$Neutralizing_antibodies
  
  # create name for model
  modtype <- gsub('R0_|_new', '', val_type)
  modtype <- paste0(toupper(substr(modtype, 1, 1)), substr(modtype, 2, nchar(modtype)), ' model')
  
  # summarize data
  sero <- sero_full %>%
    # filter(anc < 0.2) %>%
    group_by(Country, NeutAnti) %>%
    summarise(S_median = median(Seroprevalence)
              , S_5th = quantile(Seroprevalence, 0.25)
              , S_95th = quantile(Seroprevalence, 0.75)
              , R0_median = median(median)
              , R0_5th = quantile(median, 0.25)
              , R0_95th = quantile(median, 0.75)) %>%
    mutate('Model' = modtype) %>%
    as.data.frame()
  return(sero)
}

# library(viridis)
validation_plot <- function(df){
  p <- df %>%
    ggplot(aes(x = S_median, y = R0_median, group = Model)) +
    geom_errorbarh(aes(xmin = S_5th, xmax = S_95th, col = Country)) +
    geom_errorbar(aes(ymax = R0_95th, ymin = R0_5th, col = Country)) +
    geom_point(aes(x = S_median, y = R0_median)) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', fullrange = TRUE) +
    stat_cor(aes(label = after_stat(rr.label)), vjust = 1) +
    facet_grid(~ Model) +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    xlab('Seroprevalence') +
    ylab(expression(paste('Model estimated  ', R[0]))) +
    ylim(0,2) +
    theme(legend.position = 'bottom') #+
    guides(color = 'none')
  
  return(p)
  
}

cdf <- create_summary_dataset(val_type = 'R0_climate_new')
adf <- create_summary_dataset(val_type = 'R0_ancestry_new')
fdf <- create_summary_dataset(val_type = 'R0_full_new')

val.df <- do.call(rbind, list(cdf, adf, fdf))
valAll <- validation_plot(df = val.df)
ggsave(valAll, filename = 'figures/validation_all_serosurveys.pdf', width = 10, height = 5)

val.df.neutr <- subset(val.df, NeutAnti == 'Yes')
valNeut <- validation_plot(df = val.df.neutr) 
ggsave(valNeut, filename = 'figures/validation_neutralizing_serosurveys.pdf', width = 10, height = 5)
# ggsave(valNeut, filename = 'figures/validation_neutralizing_serosurveys.eps', width = 10, height = 5)

odf <- create_summary_dataset(val_type = 'R0_ancestry_omega_new')
pdf <- create_summary_dataset(val_type = 'R0_ancestry_pMI_new')
a2.df <- rbind(odf, pdf, adf)
a2.df$Model[a2.df$Model == 'Ancestry_omega model'] <- 'Omega component model' 
a2.df$Model[a2.df$Model == 'Ancestry_pMI model'] <- 'pMI component model' 

a2plot <- validation_plot(df = a2.df)
ggsave(a2plot, filename = 'figures/validation_all_serosurveys_ancestry_components.pdf', width = 10, height = 5)

a3.df <- subset(a2.df, NeutAnti == 'Yes')
a3plot <- validation_plot(df = a3.df)
ggsave(a3plot, filename = 'figures/validation_neutralizing_serosurveys_ancestry_components.pdf', width = 10, height = 5)

x <- MPP %>%
  group_by(Country, Neutralizing_antibodies, Physicians_per_100000, Childhood_immunization_coverage) %>%
  summarise(S_median = median(Seroprevalence)
            , S_5th = unname(quantile(Seroprevalence, 0.25))
            , S_95th = unname(quantile(Seroprevalence, 0.75))
            ) %>%
  as.data.frame()

s1 <- ggplot(MPP) +
  geom_point(data = x, aes(x = S_median, y = Physicians_per_100000, col = Neutralizing_antibodies), size = 3) +
  scale_color_manual(values = c('Yes' = 'darkgreen', 'No' = 'green')) +
  # geom_errorbarh(aes(xmin = 'S_5th', xmax = 'S_95th')) +
  theme_bw() +
  ylab('Physicians/100,000 population') +
  xlab('Seroprevalence (median)') +
  # theme(legend.position = "none") +  # Remove default legend
  theme(legend.position = c(0.25,0.85)) +
  guides(color = guide_legend(title = "Neutralizing antibodies"))

s2 <- ggplot(MPP) +
  geom_point(data = x, aes(x = S_median, y = Childhood_immunization_coverage, col = Neutralizing_antibodies), size = 3) +
  scale_color_manual(values = c('Yes' = 'darkgreen', 'No' = 'green')) +
  theme_bw() +
  ylab('Childhood immunizations coverage') +
  xlab('Seroprevalence (median)') +
  guides(color = "none")

healthsystem <- ggarrange(s1, s2, ncol = 2)
ggsave(filename = 'figures/Health_system_vs_seroprevalence.pdf', plot = healthsystem, width = 10, height = 5)
