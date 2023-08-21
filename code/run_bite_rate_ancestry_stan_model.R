# load libraries
library(rstan)
library(matrixStats)
library(tidyverse)

# generic traits, currently data generated from chatGPT
alpha_vals <- c(2.4, 1.7, 2.5, 1.6, 1.6, 1.7, 2.8, 2.7, 3.4, 1.9, 1.6) # skipping extreme value from BG traps (all others are HLC) 
b_vals <- c(0.80, 0.23, 0.10, 0.50, 0.33, 0.13, 0.52, 0.28, 0.23, 0.51)
EIR_vals <- c(8.7, 9.2, 11.7, 7.8, 8.8, 10, 9.3, 7.5, 8.8, 7.9)
lf_vals <- c(19.2, 17.8, 16.4, 12.1, 16.4, 17.2)
pMI_vals <- c(0.89, 0.60, 0.40, 0.82, 0.75, 0.46, 0.77, 0.55, 0.60, 0.86)

traits_prior <- list(alpha_vals, b_vals, EIR_vals, lf_vals, pMI_vals)
lapply(traits_prior, mean)
lapply(traits_prior, sd)

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
    # , ancestry_N_new = length(seq(0, 1, 0.05))
    # , aa_ancestry_new = seq(0, 1, 0.05)
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
  if(param_name != 'R0_ancestry'){
    param_name_new <- paste0(param_name, '_ancestry_new')
  } else {
    param_name_new <- param_name
  }

  samps <- data.frame(list_of_draws[param_name_new])
  sampMeans <- colMeans(samps, na.rm = T)
  sampQuantiles <- colQuantiles(list_of_draws[[param_name_new]], na.rm = T, probs = c(0.025, 0.975))
  sampQuantiles <- ifelse(sampQuantiles < 0, 0, sampQuantiles)
  
  yMin = 0
  yMax = max(sampQuantiles)
  plot(df[['aa_ancestry_new']], sampMeans, type='l', lwd=2, ylab=param_name, xlab='Proportion aa ancestry', ylim = c(0,yMax), xlim = c(0,1))
  lines(df[['aa_ancestry_new']], sampQuantiles[,1], lty=2, col='red', ylim=c(0,yMax))
  lines(df[['aa_ancestry_new']], sampQuantiles[,2], lty=2, col='red', ylim=c(0,yMax))
  if(param_name != 'R0_ancestry'){
    points(df[['aa_ancestry']], df[['omega_ancestry']], pch = 16)
  }
}

pdf('figures/bite_rate_given_ancestry_fit.pdf', width = 11, height = 8.5)
plotSamples(mod = stan_model_fit_ancestry_omega, param_name = 'omega', df = model_data)
dev.off()

pdf('figures/bite_rate_given_ancestry_R0.pdf', width = 11, height = 8.5)
plotSamples(mod = stan_model_fit_ancestry_omega, param_name = 'R0_ancestry', df = model_data)
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
