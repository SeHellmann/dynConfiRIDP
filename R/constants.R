# "2DSDT", "DDMConf" ??
DYNWEV_MODELS <- c("2DSD", "dynWEV", "dynaViTE")
RM_MODELS <- c("IRM", "PCRM", "IRMt", "PCRMt")
MODELS <- c(DYNWEV_MODELS, RM_MODELS)
DYNWEV_PARNAMES <- c("a", "v", "t0", "d", "sz", "sv", "st0", "z", "tau", "lambda", "w", "muvis", "sigvis", "svis", "s")
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

#' @keywords internal
#' @noRd
fit_context <- function(args) {
  list(
    #### initial inputs
    data = args$data,
    model = args$model,
    optim_method = args$optim_method,
    fixed = args$fixed,
    manipulations = args$manipulations,
    n_ratings = args$n_ratings,
    restr_tau = args$restr_tau,
    sym_thetas = args$sym_thetas,
    precision = args$precision,
    opts = args$opts,
    grid_search = args$grid_search,
    logging = args$logging,
    parallel = args$parallel,
    n_cores = args$n_cores,
    #### placeholders for derived values
    #### setup_logging
    log_file = NULL,
    data_file = NULL,
    #### get_model_params
    model_type = NULL,
    parnames = NULL,
    fixed_parnames = NULL,
    manipulated_parnames = NULL,
    const_parnames = NULL,
    simult_conf = NULL,
    #### process_input_data
    dependent_vars = NULL,
    initial_n_ratings = NULL,
    used_ratings = NULL,
    maxt0 = NULL,
    thetas_parnames = NULL,
    #### build_model_matrix
    model_matrix = NULL,
    beta_map = NULL,
    beta_names = NULL
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
