# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Example data
set.seed(4523)

data <- newdf[,c('Site', 'treatment', 'bio8_20', 'aaa2015')]
# data <- newdf[,c('Site', 'treatment', 'Lat', 'Lon')]
colnames(data) <- c('Site', 'treat', 'covariate1', 'covariate2')
# data <- subset(data, Site != 'Lambarene' | Site != 'Lare' | Site != 'Oyo State')
# data$covariate1 <- scale(data$covariate1)
data$covariate2 <- log(data$covariate2 + 0.0001)
data <- data %>% 
  mutate_at(c('covariate1', 'covariate2'), ~(scale(.) %>% as.vector))

data$id <- seq(1,nrow(data),1)
# Separate treated and control units
treated <- data %>% filter(treat == 1)
control <- data %>% filter(treat == 0)

# Function to calculate Euclidean distance
euclidean_dist <- function(treated_row, control_row) {
  sqrt((treated_row$covariate1 - control_row$covariate1)^2 +
         (treated_row$covariate2 - control_row$covariate2)^2)
}

# Find nearest control unit for each treated unit
matched_pairs <- treated %>%
  rowwise() %>%
  mutate(
    match_id = control$id[which.min(sapply(1:nrow(control), function(i) euclidean_dist(cur_data(), control[i, ])))],
    distance = min(sapply(1:nrow(control), function(i) euclidean_dist(cur_data(), control[i, ])))
  ) %>%
  left_join(control, by = c("match_id" = "id"), suffix = c("_treated", "_control"))

# Create the plot
ggplot() +
  geom_point(data = matched_pairs, aes(x = covariate1_treated, y = covariate2_treated), color = "blue", size = 3) +
  geom_point(data = matched_pairs, aes(x = covariate1_control, y = covariate2_control), color = "red", size = 3) +
  geom_segment(data = matched_pairs, aes(x = covariate1_treated, y = covariate2_treated, 
                                         xend = covariate1_control, yend = covariate2_control), color = "grey") +
  labs(title = "Matched Pairs Based on Euclidean Distance", x = "Covariate 1", y = "Covariate 2") +
  theme_minimal() +
  theme(legend.position = "none")

statmatch <- matched_pairs[,c('Site_treated', 'Site_control')]

statmatch_new <- statmatch %>%
  rename('Site' = 'Site_control') 

# MPP <- read.csv('../VBD-data/seroSites_Aaa.csv')
# MPP$Match_Site <- paste(MPP$Site, MPP$Country, sep = ', ')
# r0_modeled <- read.csv('../VBD-data/New_R0_values.csv') #  '../VBD-data/R0_all_models_sites.csv'

statmatch_new <- merge(statmatch_new, MPP[,c('Site', 'Lat', 'Lon', 'bio8_20', 'aaa2015', 'Seroprevalence', 'Neutralizing_antibodies')], by = 'Site')

statmatch_new <- statmatch_new %>%
  rename('Match_Site' = 'Site', 'site' = 'Site_treated')

x <- r0_modeled %>%
  left_join(statmatch_new, by = 'site')

x$Seroprevalence[x$site=='CapeVerde'] <- 10.9
x$Seroprevalence[x$site=='Ouaga'] <- 27.72

plot(x$Seroprevalence, x$Full_median, pch = 16)
plot(x$Seroprevalence[x$Neutralizing_antibodies=='Yes'], x$Full_median[x$Neutralizing_antibodies=='Yes'], pch = 16)

x2 <- subset(x, Neutralizing_antibodies=='Yes')

library(ggrepel)
library(ggpubr)

plotWithUncertainty <- function(df, xval, yval){
  yupper <- gsub('_median', '_upper', yval)
  ylower <- gsub('_median', '_lower', yval)
  
  # Calculate regression line and R-squared value
  fit <- lm(df[, yval] ~ df[, xval])
  rsquared <- summary(fit)$r.squared
  
  ggplot(df, aes_string(x = xval, y = yval)) +
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


p4 <- plotWithUncertainty(df = x, xval = 'Seroprevalence', yval = 'Climate_median')
p5 <- plotWithUncertainty(df = x, xval = 'Seroprevalence', yval = 'Ancestry_median')
p6 <- plotWithUncertainty(df = x, xval = 'Seroprevalence', yval = 'Full_median')
# gplots <-
ggarrange(p4, p5, p6, ncol = 3, nrow = 1)
