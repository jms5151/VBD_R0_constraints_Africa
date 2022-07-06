# function to predict R0 with modeled pMI for specific locations across ancestry gradient
pred_by_aaa <- function(location_code){
  # create sequence of ancestry values
  aaa <- seq(from = 0, to = 1, by = 0.01)
  # get temperature of city 
  temp_loc <- survey_data$bio.bio8_temp_wetq[survey_data$location_code == location_code]
  # new df
  newDF <- data.frame(
    'anc' = aaa # this is for pMI model
    # , 'prop_aaa_ancestry' = aaa # this is for biting rate model
    , 'bio.bio8_temp_wetq' = temp_loc
    , 'Virus' = 'Senegal_2011'
    , 'Dose' = 10^5
    )
  # Predict pMI_ancestry
  pMI_ancestry <- predict(object = infection_model, newdata = newDF, type = 'response')
  # pMI_ancestry <- pred_pMI_ancestry(df = newDF, location_code = location_code)
  # Predict biting rate reduction
  # biting_rate_reduction <- predict(object = biting_rate_reduction_model, newdata = newDF)
  # Predict R0 with predicted pMI
  lapply(pMI_ancestry, function(x) R0_NGM_adj_pMI(temp = temp_loc, pMI_ancestry = x))
  # Predict R0 with predicted pMI and adjusted biting rate
  # lapply(pMI_ancestry, function(x)
  #   R0_NGM_adj_pMI_biting_rate(
  #     temp = temp_loc
  #     , a_adjustment = biting_rate_reduction
  #     , pMI_ancestry = x))
}


# function to predict R0 with modeled pMI for specific location by year
source('code/model_biting_rate_given_ancestry.R')
source('code/function_predict_pMI_from_ancestry.R')

ts_pred_by_year <- function(location_code){
  # get temperature of city 
  temp_loc <- survey_data$bio.bio8_temp_wetq[survey_data$location_code == location_code]
  # get ancestry of city 
  aaa_value <- survey_data$prop_aaa_ancestry[survey_data$location_code == location_code]
  aaa_loc <- data.frame('prop_aaa_ancestry' = aaa_value)
  # Predict pMI_ancestry
  pMI_ancestry <- pred_pMI_ancestry(df = ts_ancestry, location_code = location_code)
  # Predict biting rate reduction
  biting_rate_reduction <- predict(object = biting_rate_reduction_model, newdata = aaa_loc)
  # Predict R0 with predicted pMI
  # lapply(pMI_ancestry, function(x) R0_NGM_adj_pMI(temp = temp_loc, pMI_ancestry = x))
  # Predict R0 with predicted pMI and adjusted biting rate
  lapply(pMI_ancestry, function(x)
    R0_NGM_adj_pMI_biting_rate(
      temp = temp_loc
      , a_adjustment = biting_rate_reduction
      , pMI_ancestry = x))
}

