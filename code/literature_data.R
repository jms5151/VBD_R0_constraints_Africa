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
                        , 'bio15_20km_precip_cv'
                        , 'Latitude'
                        , 'Longitude'
                        )
                    ], c('Site'
                         , 'dens20'
                         , 'aaa2015'
                         , 'bio8_20'
                         , 'bio15_20'
                         , 'Lat'
                         , 'Lon'
                         ))
df1$treatment <- 1

# add characteristic data from the literature
df2 <- read.csv('../VBD-data/seroSites_Aaa.csv')
df2$treatment <- 0
# df2 <- subset(df2, Country != 'Madagascar')
# not sure if this is a good idea!
# df2 <- subset(df2, Neutralizing_antibodies == 'Yes')

# combine data
# BIO15 = Precipitation Seasonality (Coefficient of Variation) - the one used in previous study
# BIO18 = Precipitation of Warmest Quarter
newdf <- rbind(df1, df2[,colnames(df1)])
# newdf$bio8_20_scaled <- scale(newdf$bio8_20)
# newdf$aaa2015_scaled <- scale(newdf$aaa2015)

# Define the sets of covariates
# all_covariates <- c('dens20', 'bio15_20', 'bio8_20') # , 'aaa2015'
# direct_covariates <- c('bio8_20_scaled', 'aaa2015_scaled')
direct_covariates <- c('bio8_20', 'aaa2015')
geog_covariates <- c('Lat', 'Lon')

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

rN = 3
# check_stat_match(cNames = all_covariates, d = 'glm', N = rN) 
# check_stat_match(cNames = all_covariates, d = 'mahalanobis', N = rN)

check_stat_match(cNames = direct_covariates, d = 'euclidean', N = rN)
# check_stat_match(cNames = direct_covariates, d = 'mahalanobis', N = rN)

# check_stat_match(cNames = geog_covariates, d = 'glm', N = rN)
# check_stat_match(cNames = geog_covariates, d = 'mahalanobis', N = rN)
# check_stat_match(cNames = geog_covariates, d = 'euclidean', N = rN)

# Function to create statistical matches
create_matches <- function(covariates, d, N) {
  # Create matches
  match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(covariates, collapse=" + "))),
                       data = newdf
                       , method = 'nearest', distance = d
                       , ratio = N
                       , replace = TRUE
                       )
  
  # Return matched data
  match_data <- match_obj$match.matrix
  matchedData <- as.data.frame(matrix(NA, nrow = nrow(match_data), ncol = N))
  for (i in 1:N) {
    matchedData[, i] <- paste0(df2[as.numeric(match_data[, i]), 'Site'], ', ', df2[as.numeric(match_data[, i]), 'Country'])
  }
  colnames(matchedData) <- c('FirstMatch', 'SecondMatch', 'ThirdMatch')
  matchedData <- data.frame('site' = df1$Site, matchedData)
  matchedData[,2:(N+1)] <- lapply(matchedData[,2:(N+1)], function(x) gsub('NA, NA', NA, x))
  return(matchedData)
}

# Function to plot diagnostic plots
# plot_diagnostics <- function(covariates, d, N) {
#   match_obj <- matchit(as.formula(paste0('treatment ~ ', paste(covariates, collapse=" + "))),
#                        data = newdf, method = 'nearest', distance = d,
#                        ratio = N, replace = TRUE)
#   # Plot diagnostics
#   plot(match_obj, type = "qq", which.xs = as.formula(paste("~", paste(covariates, collapse=" + "))))
#   plot(match_obj, type = "ecdf", which.xs = as.formula(paste("~", paste(covariates, collapse=" + "))))
#   plot(match_obj, type = "density", which.xs = as.formula(paste("~", paste(covariates, collapse=" + "))))
# }
# 
# Create matches for all covariates
# m1 <- create_matches(all_covariates, 'mahalanobis', N = rN)

# Create matches for direct covariates
g1 <- create_matches(direct_covariates, 'glm', N = rN)
# g2 <- create_matches(direct_covariates, 'mahalanobis', N = rN)

# e1 <- create_matches(geog_covariates, 'euclidean', N = rN)
# g4 <- create_matches(geog_covariates, 'mahalanobis', N = rN)

# Plot diagnostics for all covariates
# plot_diagnostics(covariates = all_covariates, d = 'mahalanobis')

# Plot diagnostics for direct covariates
# plot_diagnostics(covariates = direct_covariates, d = 'glm', N = rN)

# intermediate step
# saveUniqueSites <- function(df, typeName){
#   x <- c(df$FirstMatch[!is.na(df$FirstMatch)], df$SecondMatch[!is.na(df$SecondMatch)], df$ThirdMatch[!is.na(df$ThirdMatch)])
#   x <- unique(x)
#   write.csv(x, paste0('../VBD-data/matched_pairs_unique_list_', typeName, '.csv'), row.names = F)
# }

# saveUniqueSites(df = m1, typeName = 'mahalanobis')
# saveUniqueSites(df = g1, typeName = 'glm')

# validation

# Function to calculate standard error for each row
calculate_standard_error <- function(row_values) {
  # Remove NA values
  row_values <- na.omit(row_values)
  # Calculate standard deviation
  sd_value <- sd(row_values)
  # Calculate number of samples
  n <- length(row_values)
  if(n == 3){
    # Calculate standard error
    se <- sd_value / sqrt(n)
  } else {
    se <- 0
  }
  return(se)
}

# read in modeled results
r0_modeled <- read.csv('../VBD-data/New_R0_values.csv') #  '../VBD-data/R0_all_models_sites.csv'
# r0_modeled <- subset(r0_modeled, site != 'Thailand')
MPP <- read.csv('../VBD-data/seroSites_Aaa.csv')
MPP$Match_Site <- paste(MPP$Site, MPP$Country, sep = ', ')

# maybe boxplot of countries and then R0 points overlaid?
# test <- MPP %>%
#   group_by(Country) %>%
#   summarise(mean_seroprevalence = mean(Seroprevalence), sd_seroprevalence = sd(Seroprevalence)) %>%
#   as.data.frame()
# 
# boxplot(MPP$Seroprevalence~MPP$Country)


create_validation_df <- function(statmatch){
  
  # replace site names with seroprevalence values
  for(i in 1:nrow(MPP)){
    statmatch <- data.frame(lapply(statmatch, function(x) {
      gsub(MPP$Match_Site[i], MPP$Seroprevalence[i], x)
    }))
  }
  
  statmatch[,2:ncol(statmatch)] <- lapply(statmatch[,2:ncol(statmatch)],  as.numeric)
  
  # join datasets
  vr0 <- r0_modeled %>%
    left_join(statmatch)
  
  # fill in sites with seroprevalence from site locations
  vr0[vr0$site == 'Ouaga',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 27.72 # placeholder 
  vr0[vr0$site == 'CapeVerde',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 10.9
  # vr0[vr0$site == 'Thailand',c('FirstMatch', 'SecondMatch', 'ThirdMatch')] <- 26.57 # average of neutralizing#26.27#23.5# ##average of all studies # placeholder
  
  # calculate means and standard errors for matches for each site
  vr0$meanPrev <- rowMeans(vr0[,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], na.rm = T)
  vr0$medPrev <- apply(vr0[, c('FirstMatch', 'SecondMatch', 'ThirdMatch')], 1, median, na.rm = TRUE)
  # Apply the function to each row
  vr0$sePrev <- apply(vr0[,c('FirstMatch', 'SecondMatch', 'ThirdMatch')], 1, calculate_standard_error)
  vr0$sePrev_Upper <- vr0$meanPrev + vr0$sePrev
  vr0$sePrev_Lower <- vr0$meanPrev - vr0$sePrev
 
  return(vr0) 
}

# m1val <- create_validation_df(statmatch = m1)
g1val <- create_validation_df(statmatch = g1)
# e1val <- create_validation_df(statmatch = e1)

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
    # geom_errorbarh(aes(xmax = df[, xupper], xmin = df[, xlower]), col = 'darkgrey') +
    geom_errorbar(aes(ymax = df[, yupper], ymin = df[, ylower]), col = 'darkgrey') +
    geom_point(size = 2, color = 'black') +
    geom_smooth(method = 'lm', se = FALSE, color = 'black', fullrange = TRUE) + # Extend the line
    annotate("text", x = 2, y = 3.4, # Move text to upper left
             label = bquote(paste("R"^2, " = ", .(sprintf("%.2f", rsquared)))), 
             hjust = 0, vjust = 1, col = 'black') +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    xlab('Seroprevalence') +
    ylab(expression(paste('Model estimated  ', R[0]))) +
    ylim(0,3.5) +
    geom_text_repel(aes(label = site)) +
    theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm")) +
    ggtitle(gsub('_median', ' model', yval))
}

matchType = 'ThirdMatch'
# m1val <- subset(m1val, site != 'Thailand')
# m1val <- subset(m1val, site != 'CapeVerde')
# p1 <- plotWithUncertainty(df = m1val, xval = matchType, yval = 'Climate_median')
# p2 <- plotWithUncertainty(df = m1val, xval = matchType, yval = 'Ancestry_median')
# p3 <- plotWithUncertainty(df = m1val, xval = matchType, yval = 'Full_median')
# ggarrange(p1, p2, p3, ncol = 3, nrow = 1)

# g1val <- subset(g1val, site != 'CapeVerde')
# g1val <- subset(g1val, site != 'Thailand')
# g1val <- subset(g1val, site != 'Ouaga')
# g1val$new <- ifelse(is.na(g1val$FirstMatch), g1val$SecondMatch, g1val$FirstMatch)
# g1val$new <- ifelse(is.na(g1val$FirstMatch), g1val$ThirdMatch, g1val$FirstMatch)
p4 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Climate_median')
p5 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Ancestry_median')
p6 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'Full_median')
# gplots <-
ggarrange(p4, p5, p6, ncol = 3, nrow = 1)
ggsave(gplots, filename = 'figures/vaidation_all_serosurveys.pdf', width = 10, height = 5)
sum(is.na(g1$FirstMatch))
sum(is.na(g1$SecondMatch))
sum(is.na(g1$ThirdMatch))
length(unique(g1$FirstMatch))
length(unique(g1$SecondMatch))
length(unique(g1$ThirdMatch))


g3val <- subset(g3val, site != 'Thailand')
# g3val <- subset(g3val, site != 'Ouaga')
g3val <- subset(g3val, site != 'CapeVerde')
p7 <- plotWithUncertainty(df = g3val, xval = matchType, yval = 'Climate_median')
p8 <- plotWithUncertainty(df = g3val, xval = matchType, yval = 'Ancestry_median')
p9 <- plotWithUncertainty(df = g3val, xval = matchType, yval = 'Full_median')
# gplots <- 
ggarrange(p7, p8, p9, ncol = 3, nrow = 1)


compplots <- ggarrange(
  p1, p2, p3, p4, p5, p6, p7, p8, p9
  , ncol = 3
  , nrow = 3
  , labels = c('dens + temp + rain', '', '', 'temp + aaa', '', '', 'geog', '', ''))
ggsave(filename = 'figures/statistical_match_validation_all_surveys_mean_prev_no_thailand.pdf', plot = compplots, width = 10, height = 9)

p7 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'pMI_median')
p8 <- plotWithUncertainty(df = g1val, xval = matchType, yval = 'omega_median')
ggarrange(p7, p8, ncol = 2)
ggsave(filename = 'figures/statistical_match_temp_aaa_vaidation_pMI_omega.pdf', width = 10, height =  6)


g2val <- g1val[!is.na(g1val),]
matchType = 'ThirdMatch'
p4 <- plotWithUncertainty(df = g2val, xval = matchType, yval = 'Climate_median')
p5 <- plotWithUncertainty(df = g2val, xval = matchType, yval = 'Ancestry_median')
p6 <- plotWithUncertainty(df = g2val, xval = matchType, yval = 'Full_median')
ggarrange(p4, p5, p6, ncol = 3, nrow = 1)


create_match_plot <- function(matchName){
  mdf <- newdf
  mdf$id <- seq(1,nrow(mdf),1)
  mdf$match_id <- NA
  
  for(i in 1:nrow(g1)){
    if(!is.na(g1[i, matchName])){
      match_id <- which(mdf$Site == gsub(',.+', '', g1[i, matchName]))
      mdf$match_id[i] <- match_id
      all_match_ids <- which(g1[,matchName] == g1[i, matchName])
      mdf$match_id[match_id] <- list(all_match_ids) 
    }
  }
  
  mdf$treat <- mdf$treatment
  mdf$covariate1 <- mdf$bio8_20
  mdf$covariate2 <- mdf$aaa2015
  
  matched_data <- mdf[,c('id', 'Site', 'treat', 'covariate1', 'covariate2', 'match_id')]
  
  # Unnest the match_id list column to create a dataframe with one row per match
  expanded_matches <- matched_data %>%
    filter(treat == 1) %>%
    unnest(match_id)
  
  # Join with control units to get their covariates
  expanded_matches <- expanded_matches %>%
    left_join(matched_data %>% filter(treat == 0), by = c("match_id" = "id"), suffix = c("_treated", "_control")) %>%
    drop_na()
  
  # Create the plot
  ggplot(expanded_matches) +
    geom_point(data = newdf, aes(x = bio8_20, y = aaa2015), col = 'lightgrey') +
    geom_point(data = expanded_matches, aes(x = covariate1_treated, y = covariate2_treated), color = 'red', size = 2.5) +
    geom_point(data = expanded_matches, aes(x = covariate1_control, y = covariate2_control), color = 'black', size = 2.5) +
    geom_segment(data = expanded_matches, aes(x = covariate1_treated, y = covariate2_treated, 
                                              xend = covariate1_control, yend = covariate2_control), color = "grey") +
    labs(title = matchName, x = "Temperature", y = "Aaa") +
    theme_minimal() +
    theme(legend.position = "none") +
    # ylim(0,0.1) +
    # xlim(24,28)
    ylim(0,0.7) +
    xlim(18,32)
  
}


f_full <- create_match_plot(matchName = 'FirstMatch')
ggsave(filename = 'figures/First_match_scatter_full.pdf', plot = f_full, width = 7, height = 5)
s_full <- create_match_plot(matchName = 'SecondMatch')
ggsave(filename = 'figures/Second_match_scatter_full.pdf', plot = s_full, width = 7, height = 5)
t_full <- create_match_plot(matchName = 'ThirdMatch')
ggsave(filename = 'figures/Third_match_scatter_full.pdf', plot = t_full, width = 7, height = 5)

f_zoom <- create_match_plot(matchName = 'FirstMatch')
ggsave(filename = 'figures/First_match_scatter_zoom.pdf', plot = f_zoom, width = 7, height = 5)
s_zoom <- create_match_plot(matchName = 'SecondMatch')
ggsave(filename = 'figures/Second_match_scatter_zoom.pdf', plot = s_zoom, width = 7, height = 5)
t_zoom <- create_match_plot(matchName = 'ThirdMatch')
ggsave(filename = 'figures/Third_match_scatter_zoom.pdf', plot = t_zoom, width = 7, height = 5)

