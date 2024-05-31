# creating matching dataset

# load site characteristic data from this study
df1 <- read.csv('../VBD-data/combined_meta_colonies_fitted_clean_updated.csv')
df1 <- na.omit(df1)
df1 <- setNames(df1[, c('Location'
                        , 'dens20'
                        , 'prop_aaa_ancestry'
                        , 'bio.bio8_temp_wetq'
                        , 'bio15_20km_precip_cv'
                        , 'bio18_20km_precip_warmq'
                        , 'Latitude'
                        , 'Longitude'
                        )
                    ], c('Site'
                         , 'dens20'
                         , 'aaa2015'
                         , 'bio8_20'
                         , 'bio15_20'
                         , 'bio18_20'
                         , 'Lat'
                         , 'Lon'
                         ))
df1$treatment <- 1

# add characteristic data from the literature
df2 <- read.csv('../VBD-data/seroSites_Aaa.csv')
df2$treatment <- 0
# df2 <- subset(df2, Neutralizing_antibodies == 'Yes')
# df2 <- subset(df2, Site != 'Oyo State')
# df2 <- subset(df2, Country != 'Madagascar')

# combine data
newdf <- rbind(df1, df2[,colnames(df1)])

# add region data
regions <- read.csv('../VBD-data/regions.csv')
newdf <- merge(newdf, regions[,c('Site', 'Region')], by = 'Site')

# save
write.csv(newdf, '../VBD-data/matching_data.csv', row.names = F)

