# source model
source('code/assess_fit_infection_model_given_ancestry.R')

# load libraries
library(ggplot2)
library(ggeffects)

# Assess results and model fit
summary(infection_model)

# predict
zikv_afr_panel$y_pred <- predict(infection_model, new_data = zikv_afr_panel, type = 'response')

# observed vs predicted
cor.test(zikv_afr_panel$Infection, zikv_afr_panel$y_pred)
table(round(zikv_afr_panel$y_pred, 1), zikv_afr_panel$Infection)

# plot model fit
ggplot(zikv_afr_panel, aes(x = y_pred, y = Infection)) + 
  geom_point(alpha = .5) +
  stat_smooth(method = 'glm', se = TRUE, method.args = list(family = binomial)) +
  theme_classic() +
  xlab('Probability infection')

# plot marginal effect of ancestry
marginal_effect_ancestry <- ggpredict(infection_model, terms = 'anc')

ggplot(marginal_effect_ancestry, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  theme_classic() +
  xlab('Ancestry')

