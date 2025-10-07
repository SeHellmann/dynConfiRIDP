library(dplyr)
library(dynConfiRIDP)

data <- SATdata %>% select(-RT2, -confidence)
data <- head(data, 100)

manipulations <- list(a ~ SAT, v ~ condition)
res <- fit_rtconf_formula(
  data,
  manipulations = manipulations,
  logging = TRUE
)
print(res)
