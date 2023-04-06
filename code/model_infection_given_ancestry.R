# Model probability of mosquito infection as function of ancestry

# load libraries
library(tidyverse)
library(ggplot2)

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

# Make population and virus factored data
zikv[, c('Population', 'Virus')] <- lapply(
  zikv[, c('Population', 'Virus')]
  , function(x) as.factor(x)
)

# plots
zikv$prop_inf <- zikv$Infected/zikv$Trials
zikv$study <- c(rep('Aubry study', nrow(zikv_afr_panel)), rep('Rose study', nrow(zikv_cpv_wide)))

ggplot(zikv, aes(x = log(Dose), y = prop_inf, group = Population, name = Population, color = Population)) + 
  geom_line() +
  geom_point(size = 2) +
  facet_grid(~Virus)
  # facet_wrap(~Virus+study, nrow = 2)

# run model
# df <- subset(zikv, Virus == 'ZIKV_Senegal_2011')
# df <- subset(zikv, Virus == 'ZIKV_Cambodia_2010')
df <- zikv

mod <- glm(cbind(Infected, Trials - Infected) ~ anc + log(Dose) + Virus, family = binomial(link = 'logit'), data = df )
summary(mod)


predY <- predict(mod, newdata = df, type = 'response')

plot(df$prop_inf, predY, pch = 16, xlim = c(0,1), ylim = c(0,1), xlab = 'Observed', ylab = 'Predicted')
abline(0, 1)
