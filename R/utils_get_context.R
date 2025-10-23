#' @keywords internal
#' @noRd
get_base_context <- function(args) {
  c(args, list(
    #### placeholders for derived values
    #### setup_logging
    log_config = NULL,
    data_file_prefix = NULL,
    #### process_input_data
    dependent_vars = NULL,
    initial_n_ratings = NULL,
    used_ratings = NULL,
    maxt0 = NULL,
    #### get_model_params
    model_type = NULL,
    beta_names = NULL,
    estimated_params = NULL,
    simult_conf = NULL,
    #### build_model_matrix
    model_matrix = NULL,
    formula_params = NULL
  ))
}

#' @keywords internal
#' @noRd
get_dynwev_optimization_context <- function(context) {
  list(
    dependent_vars = context$dependent_vars,
    model_matrix = context$model_matrix,

    fixed_params = context$fixed,
    estimated_params = context$estimated_params,
    formula_params = context$formula_params,
    beta_names = context$beta_names,

    maxt0 = context$maxt0,
    restr_tau = context$restr_tau,
    precision = context$precision,
    n_ratings = context$n_ratings,
    simult_conf = context$simult_conf,
    sym_thetas = context$sym_thetas,

    optim_method = context$optim_method,
    opts = context$opts,
    logging = context$logging
  )
}