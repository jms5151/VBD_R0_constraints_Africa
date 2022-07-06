# Model probability of biting rate

# load data
survey_data <- read.csv('../data/combined_meta_colonies_fitted_clean_updated.csv')

# linear regression: this should be updated so it saturates, currently allows 
# for biting rate reduction to exceed 1
biting_rate_reduction_model <- lm(prob_biting ~ prop_aaa_ancestry, data = survey_data)

