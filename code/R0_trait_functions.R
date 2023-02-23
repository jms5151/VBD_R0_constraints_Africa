# R0 trait functions

# This is the general function for the Briere fit.
briere <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*x*(x-T0)*sqrt(Tm-x)
}

# This is the general function for the quadratic fit. 
quadratic <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*(x-T0)*(x-Tm)
}

# This is the general function for the inverted quadratic fit.
inverted_quadratic <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    24.0
  else
    1.0/(c*(x-T0)*(x-Tm))
}

# temperature dependent traits -------------------------------------------------

# Entomological parameters for the Ae. aegypti vector. 
# eggs per female per day
# EFD <- function(temp){
#   briere(temp, 8.56e-03, 14.58, 34.61)
# }

# probability egg to adult survival
# pEA <- function(temp){
#   quadratic(temp, -5.99e-03, 13.56, 38.29)
# }

# mosquito development rate (1/larval development period)
# MDR <- function(temp){
#   briere(temp, 7.86e-05, 11.36, 39.17)
# }

# biting rate
a <- function(temp){
  briere(temp, 2.02e-04, 13.35, 40.08)
}

# probability	of mosquito	infection per	bite on	an infectious	host
pMI <- function(temp){
  briere(temp, 4.91e-04, 12.22, 37.46)
}

# adult mosquito mortality rate (1/adult lifespan)
mu_m <- function(temp){
  inverted_quadratic(temp, -1.48e-01, 9.16, 37.73)
}

# parasite development rate
EIR <- function(temp){
  briere(temp, 6.65e-05, 10.68, 45.90)
}

# transmission competence: probability of human	infection	per	bite	by	an	infectious mosquito
b <- function(temp){
  briere(temp, 8.49e-04, 17.05, 35.83)
}

# Probability that a mosquito fed on DENV-infected blood becomes infected
pC <- function(temp){
  briere(temp, 4.91e-04, 12.11, 37.46)
}

# temperature independent traits -----------------------------------------------

# Human infectivity period (days) (ref: https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726)
gamma = (1/5)

# Intrinsic incubation period (days) (ref: https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726)
delta <- (1/5.9)

# Human mortality rate (Global)
# mu_h <- (1/(mean(c(71, 75))*365)) # males, females

# Human mortality rate (Africa)
mu_h <- (1/(mean(c(63, 66))*365)) # males, females