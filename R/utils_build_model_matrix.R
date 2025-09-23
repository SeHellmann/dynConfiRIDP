#' @keywords internal
#' @importFrom stats
#'   model.matrix
build_model_matrix <- function(context) {
  # add constant columns for missing manipulated parameters
  for (p in context$manipulated_parnames) {
    if (!p %in% names(context$data)) context$data[[p]] <- 1
  }

  # exit early if there are no manipulations
  if (length(context$manipulations) == 0) {
    context$model_matrix <- matrix(nrow = nrow(context$data), ncol = 0)
    context$beta_map <- setNames(integer(0), character(0))
    return(context)
  }

  individual_matrices <- list()
  all_predictor_names <- list()
  all_beta_names <- list()

  for (p in context$manipulated_parnames) {
    mat <- model.matrix(context$manipulations[[p]], context$data)
    mat_colnames <- colnames(mat)

    individual_matrices[[p]] <- mat
    all_predictor_names[[p]] <- mat_colnames
    all_beta_names[[p]] <- paste(p, mat_colnames, sep = "_")
  }

  flat_predictors <- unlist(all_predictor_names)
  flat_betas <- unlist(all_beta_names)

  is_unique <- !duplicated(flat_predictors)
  unique_predictors <- flat_predictors[is_unique]

  full_model_matrix <- do.call(cbind, individual_matrices)
  context$model_matrix <- full_model_matrix[, is_unique, drop = FALSE]

  beta_map <- match(flat_predictors, unique_predictors)
  names(beta_map) <- flat_betas
  context$beta_map <- beta_map

  context
}
