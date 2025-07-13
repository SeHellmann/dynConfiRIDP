#' @keywords internal
#' @noRd
get_model_params <- function(
  model,
  fixed = list(),
  manipulations = list()
) {
  # Define necessary parameters for the model
  if (model %in% DYNWEV_MODELS) {
    parnames <- DYNWEV_PARNAMES
    model_type <- "dynWEV"

    model_fixed_params <- switch(model,
      "2DSD" = c("lambda" = 0, "w" = 1, "sigvis" = 1, "svis" = 1),
      # "2DSDT" = c("w" = 1, "sigvis" = 1, "svis" = 1),
      "dynWEV" = c("lambda" = 0),
      "dynaViTE" = NULL
    )
  } else if (model %in% RM_MODELS) {
    parnames <- RM_PARNAMES
    model_type <- "RM"

    model_fixed_params <- if (!grepl("t", model)) {
      c("wint" = 0, "wx" = 1, "wrt" = 0)
    } else {
      NULL
    }

    model_fixed_params <- c(model_fixed_params, "rho" = if (grepl("PC", model)) -0.5 else 0)
  }

  # Combine user fixed params with defaults:
  if (!is.null(model_fixed_params)) {
    missing_fixed <- !(names(model_fixed_params) %in% names(fixed))
    fixed <- c(fixed, model_fixed_params[missing_fixed])
  }
  fixed_parnames <- names(fixed)

  # Get manipulated parameter names
  manipulated_parnames <- character(length(manipulations))

  if (length(manipulations) > 0) {
    manipulated_parnames <- vapply(manipulations, function(x) as.character(x[[2]]), character(1))

    names(manipulations) <- manipulated_parnames
  }

  const_parnames <- setdiff(parnames, c(manipulated_parnames, names(fixed)))

  # No overlaps
  if (length(intersect(const_parnames, manipulated_parnames)) > 0) {
    stop(sprintf(
      "Invalid parameter overlap: constant and manipulated parameters overlap: %s",
      paste(intersect(const_parnames, manipulated_parnames), collapse = ", ")
    ))
  }

  if (length(intersect(const_parnames, fixed_parnames)) > 0) {
    stop(sprintf(
      "Invalid parameter overlap: const and fixed parnames overlap: %s",
      paste(intersect(const_parnames, fixed_parnames), collapse = ", ")
    ))
  }

  if (length(intersect(fixed_parnames, manipulated_parnames)) > 0) {
    stop(sprintf(
      "Invalid parameter overlap: fixed and manipulated parameters overlap: %s",
      paste(intersect(fixed_parnames, manipulated_parnames), collapse = ", ")
    ))
  }

  # No coverage mismatches
  all_assigned <- c(fixed_parnames, manipulated_parnames, const_parnames)
  if (!setequal(parnames, all_assigned)) {
    missing <- setdiff(parnames, all_assigned)
    extra <- setdiff(all_assigned, parnames)
    error_msg <- "Parameter coverage mismatch:"
    if (length(missing) > 0) error_msg <- paste0(error_msg, "\n- Missing: ", paste(missing, collapse = ", "))
    if (length(extra) > 0) error_msg <- paste0(error_msg, "\n- Unexpected: ", paste(extra, collapse = ", "))
    stop(error_msg)
  }

  return(list(
    model_type = model_type,
    parnames = parnames,
    fixed = fixed,
    fixed_parnames = fixed_parnames,
    manipulations = manipulations,
    manipulated_parnames = manipulated_parnames,
    const_parnames = const_parnames
  ))
}
