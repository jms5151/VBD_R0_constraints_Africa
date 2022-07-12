source('code/model_biting_rate_given_ancestry.R')

# plot model 
pdf('figures/biting_rate_model.pdf', height = 6, width = 8)

plot(
  biting_rate_reduction_model
  , log = ''
  , pch = 16
  , ylab = 'Biting rate'
  , xlab = 'Proportion aaa ancestry'
  , main = 'Michaelis-Menten equation'
)

dev.off()

# plot predicted versus observed biting rates
# predict biting rate from model
br_data$Predicted_prob_biting <- predict(object = biting_rate_reduction_model, newdata = br_data)

# plot
pdf('figures/biting_rate_model_fit.pdf', height = 6, width = 8)

plot(
  br_data$Predicted_prob_biting
  ,   br_data$prob_biting
  , pch = 16
  , xlim = c(0,1)
  , ylim = c(0,1)
  , xlab = 'Predicted biting rate'
  , ylab = 'Observed biting rate'
)
abline(lm(br_data$Predicted_prob_biting ~ br_data$prob_biting))

dev.off()
