rm(list = ls())
pkg <- "dynConfiRIDP"

# --- CHOOSE MODE ---
# TRUE:  installs package, imports the package as a library and runs in parallel
# FALSE: uses devtools::load_all() and runs sequentially
PARALLEL <- TRUE
N_CORES <- c(2, 2)
# -------------------

# cleanup previous session
if ("future" %in% loadedNamespaces()) future::plan("sequential")
if (paste0("package:", pkg) %in% search()) detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
if (pkg %in% loadedNamespaces()) unloadNamespace(pkg)

# load package based on mode
devtools::clean_dll()
devtools::document()
if (PARALLEL) {
  devtools::install(upgrade = "never", force = TRUE)
  library(pkg, character.only = TRUE)
} else {
  devtools::load_all()
}

library(dplyr)

load("scripts/Preprocessed_SATdata.RData")
data <- SATdata %>%
  select(-RT2, -confidence) %>%
  mutate(SAT = factor(SAT), condition = factor(condition))
data <- head(data, 10000)

models <- c("dynaViTE")
fixed <- list(z = 0.5, lambda = 1, svis = 1)
manipulations <- list(a ~ SAT, v ~ condition)
n_ratings <- 6
sym_thetas <- TRUE
opts <- list(n_attempts = 3, n_restarts = 3, maxfun = 4000)

res <- fit_rtconf_models_formula(
  data,
  models = models,
  fixed = fixed,
  manipulations = manipulations,
  n_ratings = n_ratings,
  sym_thetas = sym_thetas,
  opts = opts,
  logging = TRUE,
  parallel = PARALLEL,
  n_cores = N_CORES
)

# res <- fit_rtconf_formula(
#   data,
#   model = "dynaViTE",
#   fixed = fixed,
#   manipulations = manipulations,
#   n_ratings = n_ratings,
#   sym_thetas = sym_thetas,
#   opts = opts,
#   logging = TRUE,
#   parallel = PARALLEL,
#   n_cores = N_CORES
# )

print(res)
