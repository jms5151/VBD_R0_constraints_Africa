# Model probability of mosquito infection as function of ancestry

# load library
library(tidyverse)

# load data
zikv_afr_panel <- read.delim('../VBD-data/reformatted_ZIKV_african_panel.txt')

# Make population and virus factored data
zikv_afr_panel[, c('Population', 'Virus')] <- lapply(
  zikv_afr_panel[, c('Population', 'Virus')]
  , function(x) as.factor(x)
)

# Format data, wide to long
zikv_afr_panel <- zikv_afr_panel %>%
  gather(key = 'Infection', 'N', -c(Population, Dose, Virus, anc))

zikv_afr_panel <- zikv_afr_panel[rep(seq(nrow(zikv_afr_panel)), zikv_afr_panel$N),]

# make infection status 0/1
zikv_afr_panel$Infection <- ifelse(zikv_afr_panel$Infection == 'uninfected', 0, 1)

# run model
infection_model <- glm(
  Infection ~ 0 + anc + Virus + log(Dose)
  , family = binomial(link = 'logit')
  , data = zikv_afr_panel
)
