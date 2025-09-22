#' @keywords internal
#' @importFrom stats
#'   model.matrix
build_model_matrix <- function(context) {
  # add constant columns for missing manipulated parameters
  for (p in context$manipulated_parnames) {
    if (!p %in% names(context$data)) context$data[[p]] <- 1
  }

  # prepare output containers
  model_matrix_list <- list()
  context$fit_params_cols <- list()
  context$fit_beta_parnames <- character()
  context$fit_parnames <- character()

  browser()

  for (p in context$manipulated_parnames) {
    context$model_matrix <- model.matrix(context$manipulations[[p]], context$data)
    model_matrix_colnames <- colnames(context$model_matrix)

    # save mapping of parameter to its columns
    context$fit_params_cols[[p]] <- model_matrix_colnames

    # build beta names
    context$fit_beta_parnames <- c(context$fit_beta_parnames, paste(p, model_matrix_colnames, sep = "_"))
    context$fit_parnames <- c(context$fit_parnames, rep(p, length(model_matrix_colnames)))

    model_matrix_list[[length(model_matrix_list) + 1]] <- context$model_matrix
  }

  # combine all unique columns once
  all_cols <- do.call(cbind, model_matrix_list)

  # remove duplicate columns if needed
  unique_cols <- !duplicated(colnames(all_cols))
  context$model_matrix <- all_cols[, unique_cols, drop = FALSE]

  context
}
