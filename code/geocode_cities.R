# load libraries
library(rstudioapi)
library(ggmap)

# source api key
source('../google_api_key.R')

# Load list of cities
cities <- read.csv('../VBD-data/african_cities.csv')

# format (most copied from https://en.wikipedia.org/wiki/List_of_cities_in_Africa_by_population)
cities$City <- gsub('^ ', '', cities$City)
cities$Country <- gsub('^ ', '', cities$Country)

# add lat/lon columns
cities$Longitude <- NA
cities$Latitude <- NA

# add coordinates with loop
for(i in 1:nrow(cities)){
  cityName <- paste(cities$City[i], cities$Country[i], sep = ', ')
  cityCoords <- geocode(cityName)
  cities$Longitude[i] <- cityCoords$lon
  cities$Latitude[i] <- cityCoords$lat
}

# save
write.csv('../VBD-data/african_cities_with_coordinates.csv', row.names = F)
