#' @keywords internal
#' @noRd
build_model_matrix <- function(data, manipulated_parnames, manipulations) {
  # add constant columns for missing manipulated parameters
  for (p in manipulated_parnames) {
    if (!p %in% names(data)) data[[p]] <- 1
  }

  # prepare output containers
  model_matrix_list <- list()
  fit_params_cols <- list()
  fit_beta_parnames <- character()
  fit_parnames <- character()

  for (p in manipulated_parnames) {
    model_matrix <- model.matrix(manipulations[[p]], data)
    model_matrix_colnames <- colnames(model_matrix)

    # save mapping of parameter to its columns
    fit_params_cols[[p]] <- model_matrix_colnames

    # build beta names
    fit_beta_parnames <- c(fit_beta_parnames, paste(p, model_matrix_colnames, sep = "_"))
    fit_parnames <- c(fit_parnames, rep(p, length(model_matrix_colnames)))

    model_matrix_list[[length(model_matrix_list) + 1]] <- model_matrix
  }

  # combine all unique columns once
  all_cols <- do.call(cbind, model_matrix_list)

  # remove duplicate columns if needed
  unique_cols <- !duplicated(colnames(all_cols))
  model_matrix <- all_cols[, unique_cols, drop = FALSE]

  list(
    model_matrix = model_matrix,
    fit_params_cols = fit_params_cols,
    fit_beta_parnames = fit_beta_parnames,
    fit_parnames = fit_parnames
  )
}