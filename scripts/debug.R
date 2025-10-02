library(dplyr)

devtools::load_all()

data <- SATdata %>% select(-RT2, -confidence)
data <- head(data, 100)

manipulations <- list(a ~ SAT, v ~ condition)
res <- fit_rtconf_formula(data, manipulations = manipulations, logging = TRUE)
print(res)

ratings_to_keep <- c(2, 3, 5)
data_with_gaps <- data[data$rating %in% ratings_to_keep, ]

# res <- fit_rtconf_formula(data_with_gaps, manipulations = manipulations, n_ratings = 6)
