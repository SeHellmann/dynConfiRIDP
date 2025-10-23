pkg <- "dynConfiRIDP"

if ("future" %in% loadedNamespaces()) future::plan("sequential")

if (paste0("package:", pkg) %in% search()) detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
if (pkg %in% loadedNamespaces()) unloadNamespace(pkg)

devtools::clean_dll()
devtools::document()
# devtools::load_all()

devtools::install(upgrade = "never", force = TRUE)
library(pkg, character.only = TRUE)
packageVersion(pkg)

library(dplyr)

data <- SATdata %>% select(-RT2, -confidence) %>% mutate(SAT = factor(SAT), condition = factor(condition))
models <- c("dynaViTE")
fixed <- list(z = 0.5, lambda = 1, svis = 1)
manipulations <- list(a ~ SAT, v ~ condition)
n_ratings <- 6
sym_thetas <- TRUE
opts <- list(n_attempts = 3, n_restarts = 3, maxfun = 4000)

data <- head(data, 10000)
res <- fit_rtconf_models_formula(
  data,
  models = models,
  fixed = fixed,
  manipulations = manipulations,
  n_ratings = n_ratings,
  sym_thetas = sym_thetas,
  opts = opts,
  logging = TRUE,
  parallel = TRUE,
  n_cores = c(2, 2)
)

# data <- head(data, 1000)
# res <- fit_rtconf_formula(
#   data,
#   model = "dynaViTE",
#   fixed = fixed,
#   manipulations = manipulations,
#   n_ratings = n_ratings,
#   sym_thetas = sym_thetas,
#   opts = opts,
#   logging = TRUE
# )
print(res)
