library(dplyr)

parallel <- FALSE
n_cores <- 2

if (parallel) {
  library(dynConfiRIDP)
} else {
  devtools::load_all()
}

data <- SATdata %>% select(-RT2, -confidence)
data <- head(data, 100)

manipulations <- list(a ~ SAT, v ~ condition)
res <- fit_rtconf_formula(
  data,
  manipulations = manipulations,
  logging = TRUE,
  parallel = parallel,
  n_cores = n_cores
)
print(res)
