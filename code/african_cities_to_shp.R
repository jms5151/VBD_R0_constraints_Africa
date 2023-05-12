# source function for converting points to shapefile
source('../../../R_functions/points_to_shp.R')

# load point data
cities <- read.csv('../VBD-data/african_cities_1970_2100.csv')

# format
cities <- cities[, c('City', 'Longitude', 'Latitude')]

# remove sites from Rose et al. 2015
cities <- cities[29:nrow(cities),]

# run function to convert points to shapefile and save file
pts_to_shp(df = cities
           , latName = 'Latitude'
           , lonName = 'Longitude'
           , filename = '../VBD-data/african_cities.shp')
