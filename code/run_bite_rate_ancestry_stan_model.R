# load libraries
library(rstan)
library(matrixStats)
library(tidyverse)

# load data
br_data <- read.csv('../VBD-data/combined_meta_allpops.csv')
br_data <- subset(br_data, !is.na(prop_aaa_ancestry))

# for each trait, order by temperature from low to high 
br_data <- br_data %>% arrange(prop_aaa_ancestry)

model_data <- 
  list(
    omega_ancestry_N = nrow(br_data)
    , omega_ancestry = br_data$prob
    , aa_ancestry = br_data$prop_aaa_ancestry
    , ancestry_N_new = length(seq(0, 1, 0.05))
    , aa_ancestry_new = seq(0, 1, 0.05)
  )

# fit model
stan_model_fit_ancestry_omega <- sampling(
  stan_model('code/model_biting_rate_given_ancestry.stan')
  , data = model_data
  , iter = 4000
)

# save model
saveRDS(stan_model_fit_ancestry_omega,'../models/stan_model_fit_ancestry_omega.rds')

# traceplots
pdf('figures/bite_rate_ancestry_stan_traceplots.pdf', width = 11, height = 8.5)
mod_params <- c('omega_ancestry_constant', 'omega_ancestry_d', 'omega_ancestry_e', 'omega_ancestry_sigma')
rstan::traceplot(stan_model_fit_ancestry_omega, par = c('lp__', mod_params), ncol = 2)
dev.off()

# plot model
plotSamples <- function(mod, param_name, df){
  list_of_draws <- rstan::extract(mod)
  param_name_new <- paste0(param_name, '_ancestry_new')

  samps <- data.frame(list_of_draws[param_name_new])
  sampMeans <- colMeans(samps, na.rm = T)
  sampQuantiles <- colQuantiles(list_of_draws[[param_name_new]], na.rm = T, probs = c(0.025, 0.975))
  sampQuantiles <- ifelse(sampQuantiles < 0, 0, sampQuantiles)
  
  yMin = 0
  yMax = max(sampQuantiles)
  plot(df[['aa_ancestry_new']], sampMeans, type='l', lwd=2, ylab=param_name, xlab='Proportion aa ancestry', ylim = c(0,yMax), xlim = c(0,1))
  lines(df[['aa_ancestry_new']], sampQuantiles[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(df[['aa_ancestry_new']], sampQuantiles[,2], lty=2, col='red', ylim=c(0,yMax))
  points(df[['aa_ancestry']], df[['omega_ancestry']], pch = 16)
}

pdf('figures/bite_rate_given_ancestry_fit.pdf', width = 11, height = 8.5)
plotSamples(mod = stan_model_fit_ancestry_omega, param_name = 'omega', df = model_data)
dev.off()

# paramter values
list_of_draws <- rstan::extract(stan_model_fit_ancestry_omega)
samps <- data.frame(list_of_draws$omega_ancestry_constant)
colMeans(samps, na.rm = T)

# ppc plot
# extract ppc samples and summarise
ppc_indexes <- names(stan_model_fit_ancestry_omega)[grep('ppc', names(stan_model_fit_ancestry_omega))]
ppc_estimates <- rstan::extract(stan_model_fit_ancestry_omega, ppc_indexes)
ppc_estimates <- lapply(ppc_estimates, quantile, probs=c(0.025,0.50,0.975), na.rm=TRUE)

ppc_estimates_quants <- do.call(rbind.data.frame, ppc_estimates)
colnames(ppc_estimates_quants) <- c('lower', 'median', 'upper')
ppc_estimates_quants$trait <- names(ppc_estimates)
ppc_estimates_quants$trait <- gsub('_ancestry.*', '', ppc_estimates_quants$trait)

# combine estimated and observed data
ppc_estimates_quants$value <- model_data$omega_ancestry

# plot
pdf('figures/bite_rate_ppc_plots.pdf', width = 11, height = 8.5)
ggplot(ppc_estimates_quants, aes(value, median)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_classic() + 
  facet_wrap(.~trait, scales = 'free') +
  geom_abline() +
  ggtitle('probability aa biting human host')
dev.off()