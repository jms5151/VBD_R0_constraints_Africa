# create statistical matching dataset from seroprevalence data in literature

# load libraries
library(MatchIt)
library(optmatch)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# load site characteristic data from this study
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

# add characteristic data from the literature
df2 <- read.csv('../VBD-data/seroSites_Aaa.csv')
df2$treatment <- 0

# combine data
# BIO15 = Precipitation Seasonality (Coefficient of Variation) - the one used in previous study
# BIO18 = Precipitation of Warmest Quarter
newdf <- rbind(df1, df2[,colnames(df1)])

# Define the sets of covariates
all_covariates <- c('dens20', 'bio15_20', 'bio8_20', 'aaa2015')
direct_covariates <- c('bio8_20', 'aaa2015')

# check quality of statistical matches
check_stat_match <- function(cNames, d){
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(cNames, collapse=" + "))),
                       data = newdf, method = 'nearest', distance = d,
                       ratio = 3,
                       replace = TRUE)
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

check_stat_match(cNames = all_covariates, d = 'glm')
check_stat_match(cNames = all_covariates, d = 'mahalanobis')

check_stat_match(cNames = direct_covariates, d = 'glm')
check_stat_match(cNames = direct_covariates, d = 'mahalanobis')


# Function to create statistical matches
create_matches <- function(covariates, d) {
  # Create matches
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(covariates, collapse=" + "))),
                       data = newdf, method = 'nearest', distance = d,
                       ratio = 3, replace = TRUE)
  
  # Return matched data
  match_data <- match_obj$match.matrix
  matchedData <- as.data.frame(matrix(NA, nrow = nrow(match_data), ncol = 3))
  for (i in 1:3) {
    matchedData[, i] <- paste0(df2[as.numeric(match_data[, i]), 'Site'], ', ', df2[as.numeric(match_data[, i]), 'Country'])
  }
  colnames(matchedData) <- c('FirstMatch', 'SecondMatch', 'ThirdMatch')
  matchedData <- data.frame('site' = df1$Site, matchedData)
  matchedData[,2:4] <- lapply(matchedData[,2:4], function(x) gsub('NA, NA', NA, x))
  return(matchedData)
}

# Function to plot diagnostic plots
plot_diagnostics <- function(covariates, d) {
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(covariates, collapse=" + "))),
                       data = newdf, method = 'nearest', distance = d,
                       ratio = 3, replace = TRUE)
  # Plot diagnostics
  plot(match_obj, type = "qq", which.xs = as.formula(paste("~", paste(covariates, collapse=" + "))))
  plot(match_obj, type = "ecdf", which.xs = as.formula(paste("~", paste(covariates, collapse=" + "))))
  plot(match_obj, type = "density", which.xs = as.formula(paste("~", paste(covariates, collapse=" + "))))
}

# Create matches for all covariates
m1 <- create_matches(all_covariates, 'mahalanobis')

# Create matches for direct covariates
g1 <- create_matches(direct_covariates, 'glm')

# Plot diagnostics for all covariates
plot_diagnostics(covariates = all_covariates, d = 'mahalanobis')

# Plot diagnostics for direct covariates
plot_diagnostics(covariates = direct_covariates, d = 'glm')

# intermediate step
saveUniqueSites <- function(df, typeName){
  x <- c(df$FirstMatch[!is.na(df$FirstMatch)], df$SecondMatch[!is.na(df$SecondMatch)], df$ThirdMatch[!is.na(df$ThirdMatch)])
  x <- unique(x)
  write.csv(x, paste0('../VBD-data/matched_pairs_unique_list_', typeName, '.csv'), row.names = F)
}

saveUniqueSites(df = m1, typeName = 'mahalanobis')
saveUniqueSites(df = g1, typeName = 'glm')

# validation

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

# read in modeled results
r0_modeled <- read.csv('../VBD-data/R0_all_models_sites.csv')

create_validation_df <- function(statmatch, seroSitesfilepath){
  MPP <- read.csv(seroSitesfilepath)
  
  # replace site names with seroprevalence values
  for(i in 1:nrow(MPP)){
    statmatch <- data.frame(lapply(statmatch, function(x) {
      gsub(MPP$Match_Site[i], MPP$Prevalence[i], x)
    }))
  }
  
  statmatch[,2:4] <- lapply(statmatch[,2:4],  as.numeric)
  
  # join datasets
  vr0 <- r0_modeled %>%
    left_join(statmatch)
  
  # fill in sites with seroprevalence from site locations
  vr0[vr0$site == 'Ouaga',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 27.72 # placeholder 
  vr0[vr0$site == 'CapeVerde',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 10.9
  vr0[vr0$site == 'Thailand',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 23.5#26.27 average of all studies # placeholder
  
  # calculate means and standard errors for matches for each site
  vr0$meanPrev <- rowMeans(vr0[,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], na.rm = T)
  # Apply the function to each row
  vr0$sePrev <- apply(vr0[,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], 1, calculate_standard_error)
  vr0$sePrev_Upper <- vr0$meanPrev + vr0$sePrev
  vr0$sePrev_Lower <- vr0$meanPrev - vr0$sePrev
 
  return(vr0) 
}

m1val <- create_validation_df(statmatch = m1, seroSitesfilepath = '../VBD-data/matched_pairs_prevalence_all_mahalanobis_3matches.csv')
g1val <- create_validation_df(statmatch = g1, seroSitesfilepath = '../VBD-data/matched_pairs_prevalence_glm_direct_3matches.csv')

# MPP$Prevalence <- ifelse(MPP$Match_Site=='Lare, Ethiopia', MPP$Prevalence2, MPP$Prevalence)

# thailandMeanSP <- (24.3+27.4+24.5+17+58+26.5+16.5+22.2+15.1+26.7+30.8)/11
# thailandSDSP <- sd(c(24.3,27.4,24.5,17,58,26.5,16.5,22.2,15.1,26.7,30.8))


# plotting function
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
             label = bquote(paste("R"^2, " = ", .(sprintf("%.2f", rsquared)))), 
             hjust = 0, vjust = 1, col = 'black') +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    xlab('Seroprevalence') +
    ylab(expression(paste('Model estimated  ', R[0]))) +
    geom_text_repel(aes(label = site)) +
    theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm")) +
    ggtitle(gsub('_median', ' model', yval))
}

p1 <- plotWithUncertainty(df = m1val, xval = 'meanPrev', yval = 'Climate_median')
p2 <- plotWithUncertainty(df = m1val, xval = 'meanPrev', yval = 'Ancestry_median')
p3 <- plotWithUncertainty(df = m1val, xval = 'meanPrev', yval = 'Full_median')

p4 <- plotWithUncertainty(df = g1val, xval = 'meanPrev', yval = 'Climate_median')
p5 <- plotWithUncertainty(df = g1val, xval = 'meanPrev', yval = 'Ancestry_median')
p6 <- plotWithUncertainty(df = g1val, xval = 'meanPrev', yval = 'Full_median')


compplots <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, labels = c('All covariates', '', '', 'Temperature + aaa covariates only', '', ''))

plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'omega_median')
plotWithUncertainty(df = vr0, xval = 'meanPrev', yval = 'pMI_median')
