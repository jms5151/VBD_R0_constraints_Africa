library(tidyverse)
library(rstan)

# load data
oid50 <- read.csv('../VBD-data/OID50s_ZIKV_Cambodia.csv')

# most infectious strain oid50 value 
mostInfStrain <- min(oid50$OID50[oid50$Study == 'Aubry et al.'])

# calculate pMI with reference to most infectious strain (Guadeloupe)
oid50$pMI <- mostInfStrain/oid50$OID50

# Override pMI calculation for CPV in reference to original study Guadeloupe population
oid50$pMI[oid50$Population == 'CPV' & oid50$Study == 'Rose et al.'] <- oid50$OID50[oid50$Population == 'Guadeloupe' & oid50$Study == 'Rose et al.'] / oid50$OID50[oid50$Population == 'CPV' & oid50$Study == 'Rose et al.']

# assume probability of most infectious strain is 50%, reduce all pMI values by 0.5 
oid50$pMI <- oid50$pMI - 0.5

# remove Rose et al study other than CPV
oid50 <- subset(oid50, Study == 'Aubry et al.' | Population == 'CPV')

# aaa prop
oid50$Aaa <- oid50$Aaa/100

# replace Aaa of 0 with very low value 
oid50$Aaa[oid50$Aaa == 0] <- 0.000001

# arrange by ancestry
oid50 <- oid50 %>% arrange(Aaa)

# simulate to test
sat <- function(constant, d, e, aaa){constant + ((d - constant)/(1 + (e / aaa)))}
x <- sat(0.28, 0.5, 0.05, oid50$Aaa)
plot(oid50$Aaa, oid50$pMI, pch = 16, col = 'red')
points(oid50$Aaa, x, pch = 16)

# gather data
model_data <- 
  list(
    omega_ancestry_N = nrow(oid50)
    , omega_ancestry = oid50$pMI
    , aa_ancestry = oid50$Aaa
    # , ancestry_N_new = length(seq(0, 1, 0.05))
    # , aa_ancestry_new = seq(0, 1, 0.05)
  )

# fit model
stan_model_fit_ancestry_omega <- sampling(
  stan_model('code/model_biting_rate_given_ancestry.stan')
  , data = model_data
  , iter = 2000
  # , init = list(list(omega_ancestry_constant = 0.28
  #               , omega_ancestry_d = 0.5
  #               , omega_ancestry_e = 0.05
  #               , omega_ancestry_sigma = 0.2))
  # , chains = 1
)

