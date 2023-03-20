# load libraries
library(rstan)
# library(shinystan)
# library(boot)
library(matrixStats)
library(tidyverse)

# load data
br_data <- read.csv('../VBD-data/combined_meta_allpops.csv')
br_data <- subset(br_data, !is.na(prop_aaa_ancestry))

# for each trait, order by temperature from low to high 
br_data <- br_data %>% arrange(prop_aaa_ancestry)

model_data <- 
  list(
    alpha_ancestry_N = nrow(br_data)
    , alpha_ancestry = br_data$prob
    , aa_ancestry = br_data$prop_aaa_ancestry
    , ancestry_N_new = length(seq(0, 1, 0.05))
    , aa_ancestry_new = seq(0, 1, 0.05)
  )

# fit model
stan_model_fit_ancestry_alpha <- sampling(
  stan_model('code/model_biting_rate_given_ancestry.stan')
  , data = model_data
  , iter = 4000
)

mod_params <- c('alpha_ancestry_constant', 'alpha_ancestry_d', 'alpha_ancestry_e', 'alpha_ancestry_sigma')
rstan::traceplot(stan_model_fit_ancestry_alpha, par = c('lp__', mod_params), ncol = 2)

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
  points(df[['aa_ancestry']], df[['alpha_ancestry']], pch = 16)
}

plotSamples(mod = stan_model_fit_ancestry_alpha, param_name = 'alpha', df = model_data)

list_of_draws <- rstan::extract(stan_model_fit_ancestry_alpha)
samps <- data.frame(list_of_draws$alpha_ancestry_constant)
colMeans(samps, na.rm = T)
