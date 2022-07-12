# Run alternative R0 models

# source functions and models --------------------------------------------------

# source custom functions
source('code/R0_model_functions.R')

# source component models
source('code/model_biting_rate_given_ancestry.R')
source('code/function_predict_pMI_from_ancestry.R')

# load data --------------------------------------------------------------------

# climate, mosquito ancestry, and human population data for 27 survey locations
# survey_data <- read.csv('../data/combined_meta_colonies_fitted_clean_updated.csv')
survey_data <- read.csv('../data/combined_meta_allpops.csv')
survey_data <- subset(survey_data, !is.na(prop_aaa_ancestry))

# time series data for two cities
# columns: human pop, city human pop density, ancestry (by city)
ts_ancestry <- read.delim('../data/back_envelope_aaa_over_time.txt')

# run alternative R0 models ----------------------------------------------------

# R0 model - climate only ------
R0_clim_only <- lapply(survey_data$bio.bio8_temp_wetq, function(x) R0_NGM(temp = x))

# R0 model - climate + adjusted biting rate ------
biting_rate_reduction <- predict(object = biting_rate_reduction_model, newdata = survey_data)
R0_clim_bite_rate <- R0_NGM_adj_biting_rate(temp = survey_data$bio.bio8_temp_wetq, a_adjustment = biting_rate_reduction)

# R0 model - climate + infection|ancestry ------
pMI_ancestry_pred_surveys <- pred_pMI_ancestry(df = survey_data, location_code = NA)
R0_clim_pMI <- R0_NGM_adj_pMI(temp = survey_data$bio.bio8_temp_wetq, pMI_ancestry = pMI_ancestry_pred_surveys)
# results differs considerably if prop_aaa used instead of predicted pMI_ancestry 

# R0 model - climate + infection|ancestry + transmission ------


# R0 model - climate + infection|ancestry + adjusted biting rate ------
R0_clim_pMI_bite_rate <- R0_NGM_adj_pMI_biting_rate(temp = survey_data$bio.bio8_temp_wetq
                                                    , a_adjustment = biting_rate_reduction
                                                    , pMI_ancestry = pMI_ancestry_pred_surveys)

# R0 model - climate + infection|ancestry + transmission + adjusted biting rate ------


# Plots -----

# source scatterplot function
source('code/function_R0_scatterplots.R')

# climate vs biting rate
R0_scatterplot(pred_R01 = R0_clim_only
               , pred_R02 = R0_clim_bite_rate
               , filename = 'R0_scatter_clim_vs_biting_rate'
               , xLabel = 'R0 (climate only)'
               , yLabel = 'R0 (climate & biting rate model)'
               , labelPoints = FALSE)

# climate vs pMI ancestry
R0_scatterplot(pred_R01 = R0_clim_only
               , pred_R02 = R0_clim_pMI
               , filename = 'R0_scatter_clim_vs_pMI'
               , xLabel = 'R0 (climate only)'
               , yLabel = 'R0 (climate & ancestry model)'
               , labelPoints = FALSE)

# climate vs biting rate + pMI ancestry
R0_scatterplot(pred_R01 = R0_clim_only
               , pred_R02 = R0_clim_pMI_bite_rate
               , filename = 'R0_scatter_clim_vs_pMI_biting_rate'
               , xLabel = 'R0 (climate only)'
               , yLabel = 'R0 (climate, biting rate, & ancestry model)'
               , labelPoints = TRUE)

# run models for two representative cities ------------------------------------- 

# source custom functions for location-specific predictions
source('code/functions_predict_R0_for_specific_locations.R')

# R0 vs ancestry ------ 

# predict R0 for two cities
aaa <- seq(from = 0, to = 1, by = 0.01)

R0_aaa_OGD_adj_pMI <- pred_by_aaa_adj_pMI(location_code = 'OGD')
R0_aaa_KUM_adj_pMI <- pred_by_aaa_adj_pMI(location_code = 'KUM')

R0_aaa_OGD_adj_pMI_br <- pred_by_aaa_adj_pMI_br(location_code = 'OGD')
R0_aaa_KUM_adj_pMI_br <- pred_by_aaa_adj_pMI_br(location_code = 'KUM')

# plot
pdf('figures/R0_vs_ancestry_pMI_adj_only.pdf', height = 6, width = 8)
plot(aaa, R0_aaa_OGD_adj_pMI, type = 'l', ylab = 'R0', xlab = 'Ancestry')
lines(aaa, R0_aaa_KUM_adj_pMI, type = 'l', lty = 2)
legend('topleft', bty = 'n', lty = c(1,2), legend = c('OGD', 'KUM'))
dev.off()

pdf('figures/R0_vs_ancestry_pMI_bite_rate_adj.pdf', height = 6, width = 8)
plot(aaa, R0_aaa_OGD_adj_pMI_br, type = 'l', ylab = 'R0', xlab = 'Ancestry')
lines(aaa, R0_aaa_KUM_adj_pMI_br, type = 'l', lty = 2)
legend('topleft', bty = 'n', lty = c(1,2), legend = c('OGD', 'KUM'))
dev.off()

# R0 vs year ------

# predict R0 for two cities
R0_pMI_OGD <- ts_pred_by_year(location_code = 'OGD')
R0_pMI_KUM <- ts_pred_by_year(location_code = 'KUM')

# plot results
pdf('figures/R0_vs_year_pMI_bite_rate_adj.pdf', height = 6, width = 8)
plot(ts_ancestry$Year, R0_pMI_OGD, type = 'l', ylab = 'R0', xlab = 'Year', ylim = c(0.6, 1.3))
lines(ts_ancestry$Year, R0_pMI_KUM, type = 'l', lty = 2)
legend('topleft', bty = 'n', lty = c(1,2), legend = c('OGD', 'KUM'))
dev.off()

# if only adjusting pMI
# pdf('figures/R0_vs_year_pMI_adj_only.pdf', height = 6, width = 8)
# plot(ts_ancestry$Year, R0_pMI_OGD, type = 'l', ylab = 'R0', xlab = 'Year', ylim = c(1, 2.5))
# lines(ts_ancestry$Year, R0_pMI_KUM, type = 'l', lty = 2)
# legend('topleft', bty = 'n', lty = c(1,2), legend = c('OGD', 'KUM'))
# dev.off()
