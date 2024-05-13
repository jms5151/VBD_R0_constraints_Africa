library(tidyverse)
library(rnaturalearth)

x <- read.csv('../VBD-data/R0_all_models_sites.csv')
x$Location <- x$site
xorig <- read.csv('../VBD-data/combined_meta_colonies_fitted_clean_updated.csv')

xx <- xorig[, c('Location', 'Latitude', 'Longitude')] %>%
  right_join(x[, c('Location', 'Ancestry_median', 'Climate_median', 'Full_median', 'omega_median', 'pMI_median')])

xx$Latitude[xx$Location == 'CapeVerde'] <- 16.5388
xx$Longitude[xx$Location == 'CapeVerde'] <- -23.0418

replace01 <- function(x) {
  x <- ifelse(x >= 1, 1, 0)
  x <- as.factor(x)
  return(x)
}

xx$Ancestry_median <- replace01(xx$Ancestry_median)
xx$Climate_median <- replace01(xx$Climate_median)
xx$Full_median <- replace01(xx$Full_median)

source('../google_api_key.R')

world <- ne_countries(scale='medium', returnclass = 'sf')

africa <- world %>% 
  filter(continent == "Africa")

africaBlankBig <- ggplot(data = africa) +
  geom_sf(fill = 'grey95') +
  coord_sf(xlim = c(-23, 55), ylim = c(-40, 40)) +
  theme_bw()

ggsave(africaBlankBig, file = 'figures/map_empty.jpg')

mapPoints <- function(colorVar) {
  ggplot(data = africa) +
    geom_sf(fill = 'grey95') +
    coord_sf(xlim = c(-23, 41), ylim = c(-5, 17)) +
    theme_bw() +
    geom_point(data = xx, mapping = aes_string(x = 'Longitude', y = 'Latitude', fill = colorVar, alpha = 0.5), size = 8, colour = 'black', pch = 21) +
    xlab('') +
    ylab('') +
    scale_fill_manual(
      values = c('grey30', 'red'),
      guide = guide_legend(reverse = TRUE, override.aes=list(lwd = 2))
    ) +
    theme(legend.position = 'blank') 
  
}

ancestryMap <- mapPoints(colorVar = 'Ancestry_median')
ggsave(ancestryMap, file = 'figures/map_ancestry01.jpg', width = 9.01, height = 3.5)

climateMap <- mapPoints(colorVar = 'Climate_median')
ggsave(climateMap, file = 'figures/map_climate01.jpg', width = 9.01, height = 3.5)

combinedMap <- mapPoints(colorVar = 'Full_median')
ggsave(combinedMap, file = 'figures/map_full01.jpg', width = 9.01, height = 3.5)

siteLocations <- combinedMap +
  scale_fill_manual(
    values = c('lightgrey', 'black'),
    guide = guide_legend(reverse = TRUE, override.aes=list(lwd = 2))
  ) 

ggsave(siteLocations, file = 'figures/map_sites.jpg', width = 9.01, height = 3.5)
