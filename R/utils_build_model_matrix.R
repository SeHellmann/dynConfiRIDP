#' @keywords internal
build_model_matrix <- function(context) {
  # exit early if there are no manipulations
  if (length(context$manipulations) == 0) {
    context$model_matrix <- matrix(nrow = nrow(context$data), ncol = 0)
    context$formula_params <- list()
    return(context)
  }

  # add constant columns for missing manipulated parameters
  for (p in names(context$manipulations)) {
    if (!p %in% names(context$data)) context$data[[p]] <- 1
  }

  individual_matrices <- list()
  all_predictor_names <- list()
  all_beta_names <- list()

  # 1. build individual model matrices for each manipulated parameter
  for (p in names(context$manipulations)) {
    mat <- model.matrix(context$manipulations[[p]], context$data)
    mat_colnames <- colnames(mat)

    individual_matrices[[p]] <- mat
    all_predictor_names[[p]] <- mat_colnames
    all_beta_names[[p]] <- paste(p, mat_colnames, sep = "_")
  }

  flat_predictors <- unlist(all_predictor_names)
  flat_betas <- unlist(all_beta_names)

  # 2. create one unified model_matrix with unique columns
  # this handles cases where different formulas share a predictor (e.g., v ~ condition, a ~ condition)
  # 'condition' will only appear once in the final model_matrix
  is_unique <- !duplicated(flat_predictors)
  unique_predictors <- flat_predictors[is_unique]

  full_model_matrix <- do.call(cbind, individual_matrices)
  context$model_matrix <- full_model_matrix[, is_unique, drop = FALSE]

  # 3. create the 'beta_map'.
  # this is crucial for the C++ code, it maps each beta (e.g., "v_condition")
  # to the *index* of its corresponding column in the unified model_matrix
  beta_map_indices <- match(flat_predictors, unique_predictors)
  beta_map <- split(beta_map_indices, flat_betas)

  context$beta_names <- c(context$beta_names, unique(flat_betas))

  # 4. create the final 'formula_params' list for the C++ context
  # this tells the C++ code which 'beta' indices and which 'model_matrix'
  # indices to use for calculating the value of each parameter
  context$formula_params <- list()
  for (p in names(context$manipulations)) {
    prefix <- paste0(p, "_")
    component_names <- names(beta_map)[startsWith(names(beta_map), prefix)]
    if (length(component_names) > 0) {

      context$formula_params[[p]] <- list(
        beta_indices = match(component_names, context$beta_names),
        model_matrix_indices = unlist(beta_map[component_names], use.names = FALSE)
      )
    }
  }

  context
}
