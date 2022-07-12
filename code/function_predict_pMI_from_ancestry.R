# function to predict mosquito infection from ancestry data
source('code/model_infection_given_ancestry.R')

pred_pMI_ancestry <- function(df, location_code){
  if(!is.na(location_code)){
    locationColName <- paste0(location_code, 'aaa')
    df$anc <- df[, locationColName]
  } else {
    df$anc <- df$prop_aaa_ancestry
  }
  df$Virus <- 'Senegal_2011'
  df$Dose <- 10^5
  predict(
    object = infection_model
    , newdata = df
    , type = 'response'
  )
}
