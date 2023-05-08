# source function for converting points to shapefile
source('../../../R_functions/points_to_shp.R')

# load point data
cities <- read.csv('../VBD-data/african_cities_2100aaa.csv')

# format
cities <- cities[, c('City', 'Longitude', 'Latitude')]

# run function to convert points to shapefile and save file
pts_to_shp(df = cities
           , latName = 'Latitude'
           , lonName = 'Longitude'
           , filename = '../VBD-data/african_cities.shp')
