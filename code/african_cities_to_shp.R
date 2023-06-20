# source function for converting points to shapefile
source('../../../R_functions/points_to_shp.R')

# load point data
cities <- read.csv('../VBD-data/african_cities_1970_SSP585_1970_2020_and_2090-2100.csv')

# format
cities <- cities[, c('City', 'Longitude', 'Latitude')]

# run function to convert points to shapefile and save file
pts_to_shp(df = cities
           , latName = 'Latitude'
           , lonName = 'Longitude'
           , filename = '../VBD-data/african_cities.shp')
