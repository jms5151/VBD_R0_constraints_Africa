# source R0 trait functions needed for R0 model functions below
source('code/R0_trait_functions.R')

# Base R0 model derived using next generation matrix method
R0_NGM <- function(temp){
  sqrt(
    ((a(temp) * b(temp) * PDR(temp))/(mu_m(temp)*(mu_m(temp) + PDR(temp)))) * 
      ((a(temp) * pMI(temp) * 2.0 * delta)/((delta + mu_m(temp))*(gamma + mu_h)))
  )
}

# Adjust mosquito biting rate
R0_NGM_adj_biting_rate <- function(temp, a_adjustment){
  sqrt(
    ((a_adjustment * a(temp) * b(temp) * PDR(temp))/(mu_m(temp)*(mu_m(temp) + PDR(temp)))) * 
      ((a(temp) * pMI(temp) * 2.0 * delta)/((delta + mu_m(temp))*(gamma + mu_h)))
  )
}

# Adjust probability of mosquito infection
R0_NGM_adj_pMI <- function(temp, pMI_ancestry){
  sqrt(
    ((a(temp) * b(temp) * PDR(temp))/(mu_m(temp)*(mu_m(temp) + PDR(temp)))) * 
      ((a(temp) * pMI_ancestry * 2.0 * delta)/((delta + mu_m(temp))*(gamma + mu_h)))
  )
}

# Adjust probability of mosquito infection & biting rate
R0_NGM_adj_pMI_biting_rate <- function(temp, a_adjustment, pMI_ancestry){
  sqrt(
    ((a_adjustment * a(temp) * b(temp) * PDR(temp))/(mu_m(temp)*(mu_m(temp) + PDR(temp)))) * 
      ((a(temp) * pMI_ancestry * 2.0 * delta)/((delta + mu_m(temp))*(gamma + mu_h)))
  )
}

# Adjust probability of mosquito infection & transmission
# R0_NGM_adj_pMI_trans <- function(temp, pMI_ancestry, transmission){
#   sqrt(
#     ((a(temp) * b(temp) * PDR(temp))/(mu_m(temp)*(mu_m(temp) + PDR(temp)))) * 
#       ((a(temp) * pMI_ancestry * transmission * 2.0 * delta)/((delta + mu_m(temp))*(gamma + mu_h)))
#   )
# }