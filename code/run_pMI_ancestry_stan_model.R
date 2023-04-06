# load libraries
library(rstan)
library(matrixStats)
library(tidyverse)

# load data
zikv_afr_panel <- read.delim('../VBD-data/reformatted_ZIKV_african_panel.txt')

# format
zikv_afr_panel$Trials <- zikv_afr_panel$uninfected + zikv_afr_panel$infected
colnames(zikv_afr_panel)[5] <- 'Infected'

# add 'ZIKV' to virus name
zikv_afr_panel$Virus <- paste0('ZIKV_', zikv_afr_panel$Virus)

# load data
zikv_cpv <- read.csv('../VBD-data/CPV-ZIKV.csv')

zikv_cpv_wide <- zikv_cpv %>%
  group_by(Virus, Population, Titer) %>%
  summarise(Infected = sum(Infection == 1), Trials = length(Infection)) %>%
  as.data.frame()

# rename titer as dose
colnames(zikv_cpv_wide)[3] <- 'Dose'

# add ancestry
zikv_cpv_wide$anc <- 0.23 # CPV
zikv_cpv_wide$anc[zikv_cpv_wide$Population == 'NGO'] <- 0.3738158
zikv_cpv_wide$anc[zikv_cpv_wide$Population == 'Gabon'] <-  0.073
zikv_cpv_wide$anc[zikv_cpv_wide$Population == 'Guadeloupe'] <- 1

# combine
zikv_colnames <- c('Virus' ,'Population', 'anc', 'Dose', 'Infected', 'Trials')

zikv <- rbind(
  zikv_afr_panel[, zikv_colnames],
  zikv_cpv_wide[, zikv_colnames]
)

# prepare the data
df <- zikv %>% arrange(Trials, Infected)
  
# create simulated new data
anc_new <- seq(0, 1, 0.1)
dose_new <- unname(quantile(df$Dose, c(0.25, 0.5, 0.75)))
virus_new <- c(0,1)

df_new <- expand.grid(anc_new, dose_new, virus_new)
colnames(df_new) <- c('anc', 'Dose', 'Virus')
df_new <- df_new %>% arrange(anc, Dose, Virus)

data_list <- list(
  pMI_ancestry_N = nrow(df),
  pMI_ancestry_Trials = df$Trials,
  pMI_ancestry_Infected = df$Infected,
  pMI_ancestry_X = model.matrix(~ 0 + scale(anc) + log(Dose), data = df),
  pMI_ancestry_Virus = ifelse(df$Virus == 'ZIKV_Senegal_2011', 1, 0),
  pMI_ancestry_N_new = nrow(df_new),
  pMI_ancestry_X_new = model.matrix(~0 + scale(anc) + log(Dose), data = df_new),
  pMI_ancestry_Virus_new = df_new$Virus
)

# fit the model
stan_model_fit_ancestry_pMI <- sampling(
  stan_model('code/model_infection_given_ancestry.stan')
  , data = data_list
  , iter = 1000
)

# assess fit
stan_model_fit_ancestry_pMI

pairs(stan_model_fit_ancestry_pMI)

pdf('figures/pMI_ancestry_stan_traceplots.pdf', width = 8, height = 6)
rstan::traceplot(stan_model_fit_ancestry_pMI, par = c('lp__', 'pMI_ancestry_b0', 'pMI_ancestry_beta[1]', 'pMI_ancestry_beta[2]', 'pMI_ancestry_gamma'), ncol = 2)
dev.off()

# PPC plots
ppc_indexes <- names(stan_model_fit_ancestry_pMI)[grep('prop', names(stan_model_fit_ancestry_pMI))]
ppc_estimates <- rstan::extract(stan_model_fit_ancestry_pMI, ppc_indexes)
ppc_estimates <- lapply(ppc_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)
ppc_estimates_quants <- do.call(rbind.data.frame, ppc_estimates)
colnames(ppc_estimates_quants) <- c('lower', 'median', 'upper')

ppc_estimates_quants$observed <- df$Infected
ppc_estimates_quants$trials <- df$Trials
ppc_estimates_quants$probInfObs <- ppc_estimates_quants$observed / ppc_estimates_quants$trials
  
# plot(ppc_estimates_quants$probInfObs, ppc_estimates_quants$median, pch = 16, xlim = c(0,1), ylim = c(0,1), xlab = 'Observed', ylab = 'Predicted', main = 'P(infection) | ancestry + dose + Virus')
# abline(0,1)

library(ggplot2)

ggplot(ppc_estimates_quants, aes(x = probInfObs, y = median)) +
  geom_point() +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  theme_classic() +
  geom_abline(slope=1, intercept=0) +
  ylim(0,1) +
  xlim(0,1) +
  xlab('Observed') +
  ylab('Predicted') +
  ggtitle('pMI | ancestry')

ggsave('figures/pMI_ancestry_ppc.pdf')


# Model predictions
pred_indexes <- names(stan_model_fit_ancestry_pMI)[grep('pred', names(stan_model_fit_ancestry_pMI))]
pred_estimates <- rstan::extract(stan_model_fit_ancestry_pMI, pred_indexes)
pred_estimates <- lapply(pred_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)
pred_estimates_quants <- do.call(rbind.data.frame, pred_estimates)
colnames(pred_estimates_quants) <- c('lower', 'median', 'upper')

df_new$lower <- pred_estimates_quants$lower
df_new$median <- pred_estimates_quants$median
df_new$upper <- pred_estimates_quants$upper

df_new$Virus <- ifelse(df_new$Virus == 1, 'ZIKV_Senegal_2011', 'ZIKV_Cambodia_2010')

ggplot(df_new, aes(x = anc, y = median, color = as.factor(Dose), fill = as.factor(Dose))) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, linetype = 'dashed') +
  facet_wrap(~Virus) +
  theme_classic() +
  xlab('aa ancestry') +
  ylab('Pr(Infection)')

ggsave('figures/pMI_ancestry_model.pdf')  
