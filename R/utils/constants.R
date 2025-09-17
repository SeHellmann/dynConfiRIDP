# "2DSDT", "DDMConf" ??
DYNWEV_MODELS <- c("2DSD", "dynWEV", "dynaViTE")
RM_MODELS <- c("IRM", "PCRM", "IRMt", "PCRMt")
MODELS <- c(DYNWEV_MODELS, RM_MODELS)
# "d", "muvis" ??
DYNWEV_PARNAMES <- c("a", "z", "sz", "v", "s", "sv", "tau", "muvis", "svis", "sigvis", "w", "lambda", "t0", "st0")
RM_PARNAMES <- c("a", "b", "mu1", "mu2", "s1", "s2", "rho", "wx", "wrt", "wint", "t0", "st0")
OPTIM_METHODS <- c("Nelder-Mead", "bobyqa") # "L-BFGS-B"
DEFAULT_OPTS <- list(
  n_attempts = 5,
  n_restarts = 5,
  maxfun = 8000,
  maxit = 2000,
  reltol = 1e-6,
  factr = 1e-10
)
PARALLEL_MODES <- c("none", "subject", "model", "both")

#' @keywords internal
#' @noRd
fit_context <- function(
  data,
  model,
  optim_method,
  fixed,
  manipulations,
  n_ratings,
  restr_tau,
  sym_thetas,
  precision,
  opts,
  grid_search,
  logging,
  parallel,
  n_cores
) {
  list(
    #### initial inputs
    data = data,
    model = model,
    optim_method = optim_method,
    fixed = fixed,
    manipulations = manipulations,
    n_ratings = n_ratings,
    restr_tau = restr_tau,
    sym_thetas = sym_thetas,
    precision = precision,
    opts = opts,
    grid_search = grid_search,
    logging = logging,
    parallel = parallel,
    n_cores = n_cores,
    #### placeholders for derived values
    #### setup_logging()
    logfile = NULL,
    #### get_model_params()
    model_type = NULL,
    parnames = NULL,
    fixed_parnames = NULL,
    manipulated_parnames = NULL,
    const_parnames = NULL,
    thetas_parnames = NULL,
    simult_conf = NULL,
    #### process_input_data()
    dependent_vars = NULL,
    initial_n_ratings = NULL,
    used_ratings = NULL,
    maxt0 = NULL,
    #### build_model_matrix()
    model_matrix = NULL,
    fit_params_cols = NULL,
    fit_beta_parnames = NULL,
    fit_parnames = NULL
  )
}

utils::globalVariables(c(
  "DYNWEV_MODELS",
  "RM_MODELS",
  "MODELS",
  "DYNWEV_PARNAMES",
  "RM_PARNAMES",
  "OPTIM_METHODS",
  "DEFAULT_OPTS",
  "PARALLEL_MODES"
))