#' @keywords internal
get_dynwev_optimization_context <- function(context) {
  estimated_map <- setNames(
    lapply(context$const_parnames, function(p) which(context$beta_names == p)),
    context$const_parnames
  )

  formula_map <- list()
  if (length(context$manipulated_parnames) > 0) {
    for (p in context$manipulated_parnames) {
      prefix <- paste0(p, "_")
      component_names <- names(context$beta_map)[startsWith(names(context$beta_map), prefix)]
      if (length(component_names) > 0) {
        formula_map[[p]] <- context$beta_map[component_names]
      }
    }
  }

  list(
    dependent_vars = context$dependent_vars,
    model_matrix = context$model_matrix,

    fixed_params = context$fixed,
    estimated_params = estimated_map,
    formula_params = formula_map,
    beta_names = context$beta_names,

    maxt0 = context$maxt0,
    restr_tau = context$restr_tau,
    precision = context$precision,
    n_ratings = context$n_ratings,
    simult_conf = context$simult_conf,
    sym_thetas = context$sym_thetas,

    optim_method = context$optim_method,
    opts = context$opts
  )
}