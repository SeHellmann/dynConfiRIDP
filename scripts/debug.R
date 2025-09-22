library(dplyr)

devtools::load_all()

SATdata <- SATdata %>% select(-RT2, -confidence)
SATdata <- head(SATdata, 100)

manipulations <- list(a ~ SAT, v ~ condition)

fit_rtconf_formula(SATdata, manipulations = manipulations)