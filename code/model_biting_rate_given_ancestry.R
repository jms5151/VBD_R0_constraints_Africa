# Model probability of biting humans
library(drc) # dose-response curves package

# load data
br_data <- read.csv('../VBD-data/combined_meta_allpops.csv')
br_data <- subset(br_data, !is.na(prop_aaa_ancestry))
colnames(br_data)[which(colnames(br_data) == 'prob')] <- 'prob_biting'

# Michaelis-Menten model (asymptotic)
biting_rate_reduction_model <- drm(prob_biting ~ prop_aaa_ancestry, data = br_data, fct = MM.3())

# plot results
# plot(biting_rate_reduction_model)
