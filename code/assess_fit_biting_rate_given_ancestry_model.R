# load model
source('code/model_biting_rate_given_ancestry.R')

# assess fit
summary(biting_rate_reduction_model)

# plot model vs data
plot(survey_data$prop_aaa_ancestry, survey_data$prob_biting, pch = 16, xlim = c(0, 0.8), ylim = c(0, 0.8))
abline(lm(survey_data$prob_biting ~ survey_data$prop_aaa_ancestry))