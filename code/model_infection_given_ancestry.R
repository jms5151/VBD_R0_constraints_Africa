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

# model output
summary(infection_model)

# assess model fit
# predict using same data to create model
zikv_afr_panel_predY <- predict(object = infection_model, newdata = zikv_afr_panel, type = 'response')

# turn probabilities into 0/1 response
zikv_afr_panel_predY2 <- ifelse(zikv_afr_panel_predY < 0.5, 0, 1)

# look at contingency table
xtabs(~zikv_afr_panel_predY2 + zikv_afr_panel$Infection)

# confusion matrix
library(caret)

confusionMatrix(
  as.factor(unname(zikv_afr_panel_predY2))
  , as.factor(zikv_afr_panel$Infection)
  , dnn = c("Prediction")
  )


## new data
zikv_cpv <- read.csv('../VBD-data/CPV-ZIKV.csv')

# collapse into treatments and outcomes
library(tidyverse)

zikv_cpv_wide <- zikv_cpv %>%
  group_by(Virus, Population, Titer) %>%
  summarise(Infected = sum(Infection == 1), Trials = length(Infection))

# add ancestry
zikv_cpv_wide$aaa <- 0.23 # CPV
zikv_cpv_wide$aaa[zikv_cpv_wide$Population == 'NGO'] <- 0.3738158
zikv_cpv_wide$aaa[zikv_cpv_wide$Population == 'Gabon'] <-  0.073
zikv_cpv_wide$aaa[zikv_cpv_wide$Population == 'Guadeloupe'] <- 1

## test runs to calculate OID50
test <- subset(zikv_cpv, Virus == 'ZIKV_Cambodia_2010' & Population == 'Guadeloupe')
testMod <- glm(Infection ~ log(Titer), family = binomial(link = 'logit'), data = test )
summary(testMod)
# oid50
-unname(coef(testMod)[1])/unname(coef(testMod)[2])

# calculates oid50, same as above
# library(MASS)
# dose.p(testMod)

## binomial regression
test <- subset(zikv_cpv_wide, Virus == 'ZIKV_Cambodia_2010' & Population == 'Guadeloupe')
testMod <- glm(Infected|Trials ~ log(Titer), family = binomial(link = 'logit'), data = test )
summary(testMod)
# oid50
-unname(coef(testMod)[1])/unname(coef(testMod)[2])

predY <- predict(testMod, newdata = test)
