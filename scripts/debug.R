library(dplyr)

parallel <- TRUE
n_cores <- 4

if (parallel) {
  library(dynConfiRIDP)
} else {
  devtools::load_all()
}

data <- SATdata %>% select(-RT2, -confidence)
# data <- head(data, 100)

manipulations <- list(a ~ SAT, v ~ condition)
opts <- list(n_attempts = 2, n_restarts = 2)
res <- fit_rtconf_formula(
  data,
  manipulations = manipulations,
  opts = opts,
  logging = TRUE,
  parallel = parallel,
  n_cores = n_cores
)
print(res)
