# Messy test code at the moment
library('readxl')

filepath <- '../../2022_Rose_Caldwell_R0/Raw_data_Aubry_et_al_Science_2020.xlsx'

# Get sheet names
sheet_names <- excel_sheets(filepath)

# Read all sheets to list
aubrey_data <- lapply(sheet_names, function(x)
{ as.data.frame(
  read_excel(filepath, sheet = x))
})

# Rename list elements
names(aubrey_data) <- sheet_names

zikv_afr_panel <- aubrey_data[[which(sheet_names == 'ZIKV_African_panel')]]

zikv_diss <- zikv_afr_panel[, c('Population', 'Virus', 'Dose')] 
zikv_diss <- unique(zikv_diss[complete.cases(zikv_diss), ])

# probability of transmission
zikv_mice <- aubrey_data[[which(sheet_names == 'ZIKV_mice')]]

# format columns
zikv_mice[, c('Transmission', 'logViremia')] <- lapply(
  zikv_mice[, c('Transmission', 'logViremia')]
  , function(x) as.numeric(x)
)

zikv_mice[, c('Genotype', 'Mouse', 'Population')] <- lapply(
  zikv_mice[, c('Genotype', 'Mouse', 'Population')]
  , function(x) as.factor(x)
)

# trans_model <- glm(
#   Transmission ~ 0 + Infection + Population + Virus + log(Dose)
#   , family = binomial(link = 'logit')
#   , data = zikv_afr_panel
# )