# load library
library(raster)

# import data
aaa <- raster('../VBD-data/prop_aaa_raster.tif')
bio8 <- raster('../VBD-data/wc2.1_10m_bio_8.tif')

# combine data
stacked_vars <- stack(list(aaa, bio8))

# crop data
e <- extent(-20, 55, 0, 25)
raster_crop <-crop(stacked_vars, e)
# plot(raster_crop)

# convert to dataframe
map_vars_df <- raster::as.data.frame(raster_crop, xy=TRUE)

# remove rows with missing values
map_vars_df <- map_vars_df[complete.cases(map_vars_df), ]

# save
write.csv(map_vars_df, '../VBD-data/map_variables.csv', row.names = F)
