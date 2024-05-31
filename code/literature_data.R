# create statistical matching dataset from seroprevalence data in literature

# load libraries
library(MatchIt)
library(optmatch)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# read in data
newdf <- read.csv('../VBD-data/matching_data.csv')

# Define covariates
direct_covariates <- c('bio8_20', 'aaa2015')
direct_covariates <- c('bio8_20', 'aaa2015', 'dens20', 'bio8_20', 'bio18_20', 'Lon') #, 'Region'

# check quality of statistical matches
check_stat_match <- function(cNames, d, N){
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(cNames, collapse=" + "))),
                       data = newdf
                       , method = 'nearest', distance = d
                       # , method = "nearest", distance = d, caliper = 0.2
                       , ratio = N
                       , replace = TRUE
                       )
  plot(summary(match_obj), main = cNames, xlim = c(0,1))
  summary(match_obj)
}


# assessments:
# https://cran.r-project.org/web/packages/MatchIt/vignettes/assessing-balance.html
# standardized mean difference (SMD) = 0.1, 0.05
# Variance ratios close to 1 indicate good balance because they imply the variances of the samples are similar
# variance ratios to be between .5 and 2.
# eCDF range from 0 to 1, with values closer to zero indicating better balance
# Report SMDs before and after matching for each covariate, any prognostically important interactions between covariates, and the prognostic score; this can be reported in a table or in a Love plot.
# Report summaries of balance for other statistics, e.g., the largest mean and maximum eCDF difference among the covariates and the largest SMD among squares, cubes, and interactions of the covariates.
# MatchIt provides tools for calculating each of these statistics so they can be reported with ease in a manuscript or report.

rN = 1
check_stat_match(cNames = direct_covariates, d = 'glm', N = rN)
# check_stat_match(cNames = direct_covariates, d = 'mahalanobis', N = rN)
# check_stat_match(cNames = direct_covariates, d = 'euclidean', N = rN)

newdf2 <- readRDS('../VBD-data/matching_data.RData')

# Function to create statistical matches
create_matches <- function(covariates, df, d, N, c) {
  # Create matches
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(covariates, collapse=" + "))),
                       data = df
                       , method = 'nearest'
                       , distance = d
                       # , ratio = N
                       , replace = TRUE
                       , caliper = c
                       )
  
  # Return matched data
  match_data <- match_obj$match.matrix
  matchedData <- as.data.frame(matrix(NA, nrow = nrow(match_data), ncol = N))
  for (i in 1:N) {
    matchedData[, i] <- paste0(df2[as.numeric(match_data[, i]), 'Site'], ', ', df2[as.numeric(match_data[, i]), 'Country'])
  }
  colnames(matchedData) <- c('FirstMatch')#, 'SecondMatch', 'ThirdMatch')
  matchedData <- data.frame('site' = df1$Site, matchedData)
  # matchedData[,2:(N+1)] <- lapply(matchedData[,2:(N+1)], function(x) gsub('NA, NA', NA, x))
  matchedData$site <- lapply(matchedData$site, function(x) gsub('NA, NA', NA, x))
  return(matchedData)
}

rN = 1

# ## try increasingly stringent calipers
# newids <- sample(nrow(newdf))
# df_boot <- newdf[newids,]
g1 <- create_matches(df = df_boot, df1 = df1, covariates = direct_covariates, d = 'glm', N = rN, c = 0.004)
g1$FirstMatch <- ifelse(g1$FirstMatch == 'NA, NA', NA, g1$FirstMatch)
# sites_matched <- unlist(g1$site[!is.na(g1$FirstMatch)])
# # sites_matched
# g1[!is.na(g1$FirstMatch), c('site', 'FirstMatch')]
# 
# new_df <- df_boot
# df1sub <- df1
# final_df <- data.frame()
# calipers = seq(0.005, 0.2, by = 0.0005)
# for(i in calipers){
#   if(sum(is.na(g1$FirstMatch))>0){
#     sites_matched <- unlist(g1$site[!is.na(g1$FirstMatch)])
#     final_df <- rbind(final_df, g1[!is.na(g1$FirstMatch), c('site', 'FirstMatch')])
#     new_df <- new_df[new_df$Site != sites_matched,]
#     df1sub <- df1sub[df1sub$Site != sites_matched,]
#     g1 <- create_matches(df = new_df, df1 = df1sub, covariates = direct_covariates, d = 'glm', N = rN, c = i)
#     g1$FirstMatch <- ifelse(g1$FirstMatch == 'NA, NA', NA, g1$FirstMatch)
#   }
# }
# 
# 
# 
# boot <- data.frame()
# 
# for(i in 1:3){
#   i
#   newids <- sample(nrow(newdf))
#   df_boot <- newdf[newids,]
#   g1 <- create_matches(df = df_boot, covariates = direct_covariates, d = 'glm', N = rN)
#   for(ii in 1:nrow(MPP)){
#     g1 <- data.frame(lapply(g1, function(x) {
#       gsub(MPP$Match_Site[ii], MPP$Seroprevalence[ii], x)
#     }))
#   }
#   g1val <- create_validation_df(statmatch = g1)
#   # boot <- rbind(boot, g1)  
#   g1val$FirstMatch <- as.numeric(g1val$FirstMatch)
#   p4 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Climate_median')
#   p5 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Ancestry_median')
#   p6 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Full_median')
#   ggarrange(p4, p5, p6, ncol = 3, nrow = 1)
#   
# }
# 
# boot$FirstMatch <- as.numeric(boot$FirstMatch)
# 
# bootstrapped_df <- boot %>%
#   group_by(site) %>%
#   summarise(q25 =  quantile(FirstMatch, 0.25, na.rm = T)
#             , median = quantile(FirstMatch, 0.5, na.rm = T)
#             , q75 =  quantile(FirstMatch, 0.75, na.rm = T)) %>%
#   as.data.frame()
# 
# bootstrapped_df[nrow(bootstrapped_df)+1,] <- c('CapeVerde', 10.9, 10.9, 10.9)
# bootstrapped_df[nrow(bootstrapped_df)+1,] <- c('Ouaga', 27.72, 27.72, 27.72)
# 
# x <- r0_modeled %>%
#   left_join(bootstrapped_df)

# newdf4 <- newdf %>%
#   arrange('site')
# newdf2 <- subset(newdf, Neutralizing_antibodies == 'Yes')
set.seed(10309)

calipers <- seq(0.001, 0.2, by = 0.001)
# test for # matched sites
# # unique matches
# probably can be done 10 times?
# cdf <- data.frame('Caliper_value' = NA, 'Sample' = NA, 'Num_matched_sites' = NA, 'Num_unique_matches' = NA)
cdf <- data.frame('Caliper_value' = NA, 'Num_matched_sites' = NA, 'Num_unique_matches' = NA)

# i = 0.009
for(i in 1:length(calipers)){
  # for(j in 1:10){
    newids <- sample(nrow(newdf))
    df_boot <- newdf[newids,]
    g1 <- create_matches(df = df_boot, covariates = direct_covariates, d = 'glm', N = rN, c = calipers[i])
    g1$FirstMatch <- ifelse(g1$FirstMatch == 'NA, NA', NA, g1$FirstMatch)
    sites_matched <- unlist(g1$site[!is.na(g1$FirstMatch)])
    unique_matches <- unique(unlist(g1$FirstMatch[!is.na(g1$FirstMatch)]))
    cdf[i,] <- c(calipers[i], length(sites_matched), length(unique_matches))
  #   cdf[i,] <- c(calipers[i], j, length(sites_matched), length(unique_matches))
  # }
}

plot(cdf$Caliper_value, cdf$Num_matched_sites, type = 'l', ylab = 'Number of matched sites', xlab = 'Caliper value')
plot(cdf$Caliper_value, cdf$Num_unique_matches, type = 'l', ylab = 'Unique matches', xlab = 'Caliper value')

cal_val_to_use <- cdf$Caliper_value[which(cdf$Num_matched_sites == max(cdf$Num_matched_sites))][1]
min_cal_val_to_use <- cdf$Caliper_value[which(cdf$Num_unique_matches == max(cdf$Num_unique_matches))[1]]

boot_df <- data.frame('Iter' = NA, 'Climate_R2' = NA, 'Ancestry_R2' = NA, 'Full_R2' = NA)

calcR2 <- function(df, yval, xval){
  fit <- lm(df[, yval] ~ df[, xval])
  rsquared <- summary(fit)$r.squared
  return(rsquared)  
}

# Maybe need to do a combination of caliper and fit
# c = 0.009
for(i in 1:50){
  newids <- sample(nrow(newdf))
  df_boot <- newdf[newids,]
  g1 <- create_matches(df = df_boot, covariates = direct_covariates, d = 'glm', N = rN, c = cal_val_to_use)
  g1val <- create_validation_df(statmatch = g1)
  # Calculate regression line and R-squared value
  Cr2 <- calcR2(g1val, 'Climate_median', 'FirstMatch')  
  Ar2 <- calcR2(g1val, 'Ancestry_median', 'FirstMatch')  
  Fr2 <- calcR2(g1val, 'Full_median', 'FirstMatch')
  boot_df[i,] <- c(i, Cr2, Ar2, Fr2)
}

boot_summary <- boot_df %>%
  summarise_all(median)

# newids <- sample(nrow(newdf))
# df_boot <- newdf[newids,]
g1 <- create_matches(df = df_boot, covariates = direct_covariates, d = 'glm', N = rN, c = cal_val_to_use)
g1val <- create_validation_df(statmatch = g1)
matchType = 'FirstMatch'
p4 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Climate_median')
p5 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Ancestry_median')
p6 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Full_median')
ggarrange(p4, p5, p6, ncol = 3, nrow = 1)
