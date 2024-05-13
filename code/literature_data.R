x <- read.csv('../VBD-data/historical_prevalence_literature.csv')
x2 <- subset(x, Age_stratified == 'No')

x2$r0 <- ifelse(!is.na(x2$Adjusted_percent_positive), x2$Adjusted_percent_positive, x2$Percent_positive)
boxplot(x2$r0~x2$Country)

x2$r0 <- 1/(100-x2$r0)
boxplot(x2$r0~x2$Country)



xx <- read.csv('C:/Users/jamie/OneDrive - Princeton University/Projects/2023_Caldwell_Rose_Lamberts_The_Lancet_Planetary_Health/VBD-data/age_stratified_data.csv')
xx$age_group_start <- as.numeric(gsub('-.*', '', xx$Age_group) )
xx$age_group_end <- as.numeric(gsub('.*-| yrs', '', xx$Age_group) )
xx$Midpoint_age <- apply(xx, 1, function(x) midpoint(x['age_group_start'], x['age_group_end']))
xx$age_prev <- xx$Midpoint_age * (xx$Prevalence/100)

library(tidyverse)

test <- xx %>%
  group_by(Country, Site, Year, Life_expectency_world_bank_stats, Neutralizing) %>%
  summarise(mean_age_infection = sum(age_prev)/sum(Prevalence/100)) %>%
  mutate(R0 = (Life_expectency_world_bank_stats/mean_age_infection) + 1) %>%
  as.data.frame()

boxplot(test$R0~test$Country)
boxplot(test$R0~test$Neutralizing)
quantile(test$R0, na.rm = T)

test2 <- x2 %>%
  left_join(test)

plot(test2$r0, test2$R0, pch = 16)


modR0 <- read.csv('C:/Users/jamie/OneDrive - Princeton University/Projects/2023_Caldwell_Rose_Lamberts_The_Lancet_Planetary_Health/VBD-data/R0_all_models_sites.csv')



# for each country pyramid calculate 1. total by age group, then 2. percentage of total by age group
# then calculate population size by multiplying site population by percentage of total within each age group

# only need to run once to calculate population sums and proportions
poppyramids <- list.files('../VBD-data/population_data/', full.names = T)

for(i in poppyramids){
  p <- read.csv(i)
  p$total_sum <- rowSums(p[,c('M', 'F')])
  totalPop <- sum(p$total_sum)
  p$total_prop <- p$total_sum/totalPop
  write.csv(p, file = i, row.names = F)
}


agestrat <- read.csv('../VBD-data/age_stratified_data.csv')

mround <- function(x,base){ base*round(x/base) }

midpoint <- function(n1,n2){ median(seq(n1,n2,1)) }

agestrat$age_group_start <- as.numeric(gsub('-.*', '', agestrat$Age_group) )
agestrat$age_group_end <- as.numeric(gsub('.*-| yrs', '', agestrat$Age_group) )

agestrat$agebin1 <- mround(agestrat$age_group_start, 5)
agestrat$agebin2 <- mround(agestrat$age_group_end, 5)

agestrat$Midpoint_age <- apply(agestrat, 1, function(x) midpoint(x['age_group_start'], x['age_group_end']))

sitepop <- read.csv('../VBD-data/total_population.csv')

library(tidyverse)

x <- agestrat %>%
  left_join(sitepop)

xx <- x[,c('Country', 'Site', 'Year', 'Age_group', 'agebin1', 'agebin2', 'Total_population')] 
write.csv(xx, file = '../VBD-data/age_stratified_data_v2.csv', row.names = F)

agestrat2 <- agestrat %>%
  left_join(xx)

agestrat2$Population_size <- agestrat2$Total_population * agestrat2$Pop_proportion


#### START HERE
library(tidyverse)
agestrat2 <- read.csv('../VBD-data/age_stratified_data_v2.csv')
agestrat2$prev <- agestrat2$Prevalence/100
agestrat2$age_prev <- agestrat2$Midpoint_age * agestrat2$prev
agestrat2$age_prev_pop <- agestrat2$Midpoint_age * agestrat2$prev * agestrat2$Population_size
agestrat2$prev_pop <- agestrat2$prev * agestrat2$Population_size

# this yields some inf values (and some NaN), need to check:)
xx <- agestrat2 %>%
  group_by(Country, Site, Year, Life_expectency_world_bank_stats, Neutralizing) %>%
  reframe(age_prev_summed = sum(age_prev),
          prev_summed = sum(prev),
          mean_prev = mean(prev),
          age_prev_pop_summed = sum(age_prev_pop),
          prev_pop_summed = sum(prev_pop),
          mean_age_infection = age_prev_summed/prev_summed,
          mean_age_infection_weighted = age_prev_pop_summed/prev_pop_summed) %>%
  # mutate_all(function(x) ifelse(is.infinite(x)|is.nan(x), 999, x)) %>%
  group_by(Country, Site, Neutralizing) %>%
  mutate(R0 = Life_expectency_world_bank_stats/mean_age_infection + 1,
         R0weighted = Life_expectency_world_bank_stats/mean_age_infection_weighted + 1) %>%
  as.data.frame()

# xx$mean_age_infection <- ifelse(xx$mean_age_infection == 999, xx$Life_expectency_world_bank_stats, xx$mean_age_infection)
# xx$mean_age_infection_weighted <- ifelse(xx$mean_age_infection_weighted == 999, xx$Life_expectency_world_bank_stats, xx$mean_age_infection_weighted)

boxplot(xx$R0~xx$Country)
boxplot(xx$R0weighted~xx$Country)
boxplot(xx$R0~xx$Neutralizing)
boxplot(xx$R0weighted~xx$Neutralizing)

hist(test$R0weighted)
## need to 

totalprev <- read.csv('../VBD-data/literature_prevalence.csv')

tp <- subset(totalprev, Country == 'Burkina Faso' | Country == 'Cape Verde' | Country == 'Gabon' | Country == 'Kenya' | Country == 'Nigeria' | Country == 'Ghana' | Country == 'Senegal' )
boxplot(tp$Prevalence~tp$Country)


df1 <- read.csv('../VBD-data/combined_meta_colonies_fitted_clean_updated.csv')
df1 <- na.omit(df1)
df1 <- setNames(df1[, c('Location'
                        , 'dens20'
                        , 'prop_aaa_ancestry'
                        , 'bio.bio8_temp_wetq'
                        , 'bio15_20km_precip_cv')
                    ], c('Site'
                         , 'dens20'
                         , 'aaa2015'
                         , 'bio8_20'
                         , 'bio15_20'
                         ))

df1$treatment <- 1
df2 <- read.csv('../VBD-data/seroSites_Aaa.csv')
df2$treatment <- 0

newdf <- rbind(df1, df2[,colnames(df1)])

# Calculate Euclidean distance
library(MatchIt)
library(optmatch)
library(tidyverse)
# BIO15 = Precipitation Seasonality (Coefficient of Variation) - the one used in previous study
# BIO18 = Precipitation of Warmest Quarter
# cols1 <- c('dens20', 'bio15_20', 'bio18_20', 'aaa2015', 'Lon', 'Lat', 'treatment') 
# cols2 <- c('dens20', 'bio15_20', 'bio18_20', 'aaa2015', 'treatment') 
# cols3 <- c('dens20', 'bio15_20', 'bio18_20', 'Lon', 'Lat', 'treatment') 
# cols4 <- c('dens20', 'bio15_20', 'bio18_20', 'treatment') 

# maybe we should first add prevalence data from sites with prevalence data 
# Burkina Faso, Cape Verde,
# cNames = cols2


check_stat_match <- function(cNames, d){
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(cNames, collapse=" + "))),
                       data = newdf, method = 'nearest', distance = d,
                       ratio = 3,
                       replace = TRUE)
  plot(summary(match_obj), main = cNames, xlim = c(0,1))
  summary(match_obj)
}

check_stat_match(cNames = c('dens20', 'bio15_20', 'bio8_20', 'aaa2015'), d = 'glm')
check_stat_match(cNames = c('bio15_20', 'bio8_20', 'aaa2015'), d = 'glm')
check_stat_match(cNames = c('dens20', 'bio15_20', 'bio8_20'), d = 'glm')
check_stat_match(cNames = c('bio15_20', 'bio8_20'), d = 'glm')
check_stat_match(cNames = c('aaa2015', 'bio8_20'), d = 'glm')

check_stat_match(cNames = c('dens20', 'bio15_20', 'bio8_20', 'aaa2015'), d = 'mahalanobis')
check_stat_match(cNames = c('bio15_20', 'bio8_20', 'aaa2015'), d = 'mahalanobis')
check_stat_match(cNames = c('dens20', 'bio15_20', 'bio8_20'), d = 'mahalanobis')

all <- c('dens20', 'bio15_20', 'bio8_20', 'aaa2015')
direct <- c('bio8_20', 'aaa2015')

check_stat_match(cNames = all, d = 'glm')
check_stat_match(cNames = all, d = 'mahalanobis')

check_stat_match(cNames = direct, d = 'glm')
check_stat_match(cNames = direct, d = 'mahalanobis')

# cNames <- colnames(newdf)[!grepl('treatment|Site|aaa2015|dens20', colnames(newdf))]

# plot(match_obj, type = "density", which.xs = ~dens20 + aaa2015 + bio8_20 + bio15_20)


create_stat_match <- function(cNames, d){
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(cNames, collapse=" + "))),
                       data = newdf, method = 'nearest', distance = d,
                       ratio = 3,
                       replace = TRUE)
  
  first_match <- as.numeric(match_obj$match.matrix[,1])
  second_match <- as.numeric(match_obj$match.matrix[,2])
  third_match <- as.numeric(match_obj$match.matrix[,3])
  matchedData <- as.data.frame(cbind(
    'site' = df1$Site
    , 'FirstMatch' = paste0(df2[first_match, 'Site'], ', ', df2[first_match, 'Country'])
  # ))
    , 'SecondMatch' = paste0(df2[second_match, 'Site'], ', ', df2[second_match, 'Country'])
    , 'ThirdMatch' = paste0(df2[third_match, 'Site'], ', ', df2[third_match, 'Country'])))
  matchedData[,2:4] <- lapply(matchedData[,2:4], function(x) gsub('NA, NA', NA, x))
  # matchedData$FirstMatch <- gsub('NA, NA', NA, matchedData$FirstMatch)
  return(matchedData)
}

g1 <- create_stat_match(cNames = direct, d = 'glm')
# m1 <- create_stat_match(cNames = c('dens20', 'bio15_20', 'bio8_20', 'aaa2015'), d = 'mahalanobis')
# m2 <- create_stat_match(cNames = c('dens20', 'bio15_20', 'bio8_20'), d = 'mahalanobis')



match_obj <- matchit(as.formula('treatment ~ bio15_20 + bio8_20'),
                     data = newdf, method = 'nearest', distance = 'glm',
                     ratio = 3,
                     replace = TRUE)

m1 <- create_stat_match(cNames = all, d = 'mahalanobis')

match_obj <- matchit(as.formula('treatment ~ bio8_20 + aaa2015'),
                     # match_obj <- matchit(as.formula('treatment ~ bio15_20 + bio8_20'),
                     data = newdf, method = 'nearest', distance = 'mahalanobis',
                     ratio = 3,
                     replace = TRUE)

plot(summary(match_obj))
plot(match_obj, type = "qq", which.xs = ~dens20 + bio15_20 + bio8_20 + aaa2015)
plot(match_obj, type = "ecdf", which.xs = ~dens20 + bio15_20 + bio8_20 + aaa2015)
plot(match_obj, type = "density", which.xs = ~dens20 + bio15_20 + bio8_20 + aaa2015)


x <- unique(m1$FirstMatch[!is.na(m1$FirstMatch)])
write.csv(x, '../VBD-data/matched_pairs_unique_list.csv', row.names = F)

# t1 <- create_stat_match(cNames = cols1, d = 'glm')
# t2 <- create_stat_match(cNames = cols2, d = 'glm')
# t3 <- create_stat_match(cNames = cols3, d = 'glm')
# t4 <- create_stat_match(cNames = cols4, d = 'glm')
# 
# m1 <- create_stat_match(cNames = cols1, d = 'mahalanobis')
# m2 <- create_stat_match(cNames = cols2, d = 'mahalanobis')
# m3 <- create_stat_match(cNames = cols3, d = 'mahalanobis')
# m4 <- create_stat_match(cNames = cols4, d = 'mahalanobis')
# 
# countNAs <- function(x){
#   cat('First match =',
#     sum(is.na(x$FirstMatch)),
#     '\nSecond match = ',
#     sum(is.na(x$SecondMatch)),
#     '\nThird match =',
#     sum(is.na(x$ThirdMatch))
#   )
# }
# 
# countNAs(x = t1)
# countNAs(x = t2)
# countNAs(x = t3)
# countNAs(x = t4)
# countNAs(x = m1)
# countNAs(x = m2)
# countNAs(x = m3)
# countNAs(x = m4)


# other assessments:
# https://cran.r-project.org/web/packages/MatchIt/vignettes/assessing-balance.html
# standardized mean difference (SMD) = 0.1, 0.05
# Variance ratios close to 1 indicate good balance because they imply the variances of the samples are similar
# variance ratios to be between .5 and 2.
# eCDF range from 0 to 1, with values closer to zero indicating better balance
# Report SMDs before and after matching for each covariate, any prognostically important interactions between covariates, and the prognostic score; this can be reported in a table or in a Love plot.
# Report summaries of balance for other statistics, e.g., the largest mean and maximum eCDF difference among the covariates and the largest SMD among squares, cubes, and interactions of the covariates.
# MatchIt provides tools for calculating each of these statistics so they can be reported with ease in a manuscript or report.

# t2 is most robust, finding most first matches and most matches overall
x <- c(g1$FirstMatch[!is.na(g1$FirstMatch)], g1$SecondMatch[!is.na(g1$SecondMatch)], g1$ThirdMatch[!is.na(g1$ThirdMatch)])
x <- c(m1$FirstMatch[!is.na(m1$FirstMatch)], m1$SecondMatch[!is.na(m1$SecondMatch)], m1$ThirdMatch[!is.na(m1$ThirdMatch)])
x <- unique(x)
write.csv(x, '../VBD-data/matched_pairs_unique_list.csv', row.names = F)

validation_df <- g1
validation_df <- m1

MPP <- read.csv('../VBD-data/matched_pairs_prevalence_glm_direct_3matches.csv')
MPP <- read.csv('../VBD-data/matched_pairs_prevalence_all_mahalanobis_3matches.csv')

MPP$Prevalence <- ifelse(MPP$Match_Site=='Lare, Ethiopia', MPP$Prevalence2, MPP$Prevalence)

for(i in 1:nrow(MPP)){
  validation_df <- data.frame(lapply(validation_df, function(x) {
    gsub(MPP$Match_Site[i], MPP$Prevalence[i], x)
  }))
}

validation_df[,2:4] <- lapply(validation_df[,2:4],  as.numeric)
# validation_df$FirstMatch <- as.numeric(validation_df$FirstMatch)
r0 <- read.csv('../VBD-data/R0_all_models_sites.csv')

vr0 <- r0 %>%
  left_join(validation_df)

vr0[vr0$site == 'Ouaga',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 27.72 # placeholder 
vr0[vr0$site == 'CapeVerde',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 10.9
vr0[vr0$site == 'Thailand',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 23.5#26.27 average of all studies # placeholder
# vr0 <- subset(vr0, site != 'Thailand')
vr0$meanPrev <- rowMeans(vr0[,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], na.rm = T)


# Function to calculate standard error for each row
calculate_standard_error <- function(row_values) {
  # Remove NA values
  row_values <- na.omit(row_values)
  # Calculate standard deviation
  sd_value <- sd(row_values)
  # Calculate number of samples
  n <- length(row_values)
  # Calculate standard error
  se <- sd_value / sqrt(n)
  return(se)
}

# Apply the function to each row
vr0$sePrev <- apply(vr0[,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], 1, calculate_standard_error)
vr0$sePrev_Upper <- vr0$meanPrev + vr0$sePrev
vr0$sePrev_Lower <- vr0$meanPrev - vr0$sePrev

# vr0$initMatch <- ifelse(is.na(vr0$FirstMatch), vr0$SecondMatch, vr0$FirstMatch)
# vr0$initMatch <- ifelse(is.na(vr0$FirstMatch), vr0$ThirdMatch, vr0$FirstMatch)

vr0$min <- as.numeric(NA)
vr0$max <- as.numeric(NA)

for(i in 1:nrow(vr0)){
  vr0$min[i] <- min(vr0[i ,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], na.rm = T)
  vr0$max[i] <- max(vr0[i ,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], na.rm = T)
}

# thailandMeanSP <- (24.3+27.4+24.5+17+58+26.5+16.5+22.2+15.1+26.7+30.8)/11
# thailandSDSP <- sd(c(24.3,27.4,24.5,17,58,26.5,16.5,22.2,15.1,26.7,30.8))

simpleplot <- function(x,y){
  model <- lm(vr0[,y] ~ vr0[,x])
  r_squared <- summary(model)$r.squared
  plot(vr0[,x], vr0[,y], pch = 16, ylab = 'Model R0', xlab = 'Seroprevalence', main = paste(y, x, sep = ', '))
  abline(model)
  legend("topleft", legend = paste("R\u00B2 =", round(r_squared, 2)), bty = "n")
}

par(mfrow=c(1,3))
simpleplot('FirstMatch', 'Ancestry_median')
simpleplot('FirstMatch', 'Full_median')
simpleplot('FirstMatch', 'Climate_median')


par(mfrow=c(3,3), mai = c(0.5,0.5,0.5, 0.5))
# par(mfrow=c(3,4), mai = c(0.5,0.5,0.5,0.5))
simpleplot('FirstMatch', 'Ancestry_median')
simpleplot('SecondMatch', 'Ancestry_median')
simpleplot('ThirdMatch', 'Ancestry_median')
# simpleplot('meanPrev', 'Ancestry_median')

simpleplot('FirstMatch', 'Full_median')
simpleplot('SecondMatch', 'Full_median')
simpleplot('ThirdMatch', 'Full_median')
# simpleplot('meanPrev', 'Full_median')

simpleplot('FirstMatch', 'Climate_median')
simpleplot('SecondMatch', 'Climate_median')
simpleplot('ThirdMatch', 'Climate_median')
# simpleplot('meanPrev', 'Climate_median')


par(mfrow=c(2,3), mai = c(0.5,0.5,0.5, 0.5))
simpleplot('meanPrev', 'Ancestry_median')
simpleplot('meanPrev', 'Full_median')
simpleplot('meanPrev', 'Climate_median')

simpleplot('initMatch', 'Ancestry_median')
simpleplot('initMatch', 'Full_median')
simpleplot('initMatch', 'Climate_median')

simpleplot('meanPrev', 'Ancestry_median')
simpleplot('meanPrev', 'omega_median')
simpleplot('meanPrev', 'pMI_median')


library(tidyverse)
rm(df)
methods_vec <- c('nearest', 'full')
d_vec <- c('glm' , "mahalanobis", "euclidean")
cNames_list <- list(cols1, cols2, cols3, cols4)
names(cNames_list) <- c('full', 'noCoords', 'noAAA', 'min')

for(i in 1:length(cNames_list)){
  cNamesI <- cNames_list[[i]]
  for(j in 1:length(methods_vec)){
    mthdJ <- methods_vec[j]
    for(k in 1:length(d_vec)){
      dK <- d_vec[k]
      x <- create_stat_match(cNames = cNamesI, mthd = mthdJ, d = dK, vType = names(cNames_list)[i])
      if(length(df) > 0){
        df <- cbind(df, x[,2])
      } else {
        df <- x
      }
    }
  }
}

library(ggplot2)
library(ggrepel)

plotWithUncertainty <- function(df, xval, yval){
  xupper <- 'sePrev_Upper'
  xlower <- 'sePrev_Lower'
  yupper <- gsub('_median', '_upper', yval)
  ylower <- gsub('_median', '_lower', yval)
  
  # Calculate regression line and R-squared value
  fit <- lm(df[, yval] ~ df[, xval])
  rsquared <- summary(fit)$r.squared
  
  ggplot(df, aes_string(x = xval, y = yval)) +
    geom_errorbarh(aes(xmax = df[, xupper], xmin = df[, xlower]), col = 'darkgrey') +
    geom_errorbar(aes(ymax = df[, yupper], ymin = df[, ylower]), col = 'darkgrey') +
    geom_point(size = 2, color = 'black') +
    geom_smooth(method = 'lm', se = FALSE, color = 'black', fullrange = TRUE) + # Extend the line
    annotate("text", x = min(df[, xval]), y = max(df[, yupper]), # Move text to upper left
             label = paste("R-squared = ", sprintf("%.2f", rsquared)), 
             hjust = 0, vjust = 1, col = 'black') +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    xlab('Seroprevalence') +
    ylab(expression(paste('Model estimated  ', R[0]))) +
    geom_text_repel(aes(label = site)) +
    theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm")) +
    ggtitle(gsub('_median', ' model', yval))
}

plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'Climate_median')
plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'Ancestry_median')
plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'Full_median')

plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'omega_median')
plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'pMI_median')
