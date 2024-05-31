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

newdf <- readRDS('../VBD-data/matching_data.RData')

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

g1 <- create_matches(df = newdf, covariates = direct_covariates, d = 'glm', N = rN, c = cal_val_to_use)
g1val <- create_validation_df(statmatch = g1)
matchType = 'FirstMatch'
p4 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Climate_median')
p5 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Ancestry_median')
p6 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Full_median')
ggarrange(p4, p5, p6, ncol = 3, nrow = 1)
