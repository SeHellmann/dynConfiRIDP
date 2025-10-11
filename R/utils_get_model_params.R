#' @keywords internal
get_model_params <- function(context) {
  # Define necessary parameters for the model
  if (context$model %in% DYNWEV_MODELS) {
    context$model_type <- "dynWEV"
    parnames <- DYNWEV_PARNAMES

    model_fixed_params <- switch(context$model,
      "2DSD" = c("lambda" = 0, "w" = 1, "sigvis" = 1, "svis" = 1),
      # "2DSDT" = c("w" = 1, "sigvis" = 1, "svis" = 1),
      "dynWEV" = c("lambda" = 0),
      "dynaViTE" = NULL
    )
  } else if (context$model %in% RM_MODELS) {
    context$model_type <- "RM"
    parnames <- RM_PARNAMES

    model_fixed_params <- if (!grepl("t", context$model)) {
      c("wint" = 0, "wx" = 1, "wrt" = 0)
    } else {
      NULL
    }

    model_fixed_params <- c(model_fixed_params, "rho" = if (grepl("PC", context$model)) -0.5 else 0)
  }

  # Combine user fixed params with defaults:
  if (!is.null(model_fixed_params)) {
    missing_fixed <- !(names(model_fixed_params) %in% names(context$fixed))
    context$fixed <- c(context$fixed, model_fixed_params[missing_fixed])
  }
  fixed_parnames <- names(context$fixed)

  # Get manipulated parameter names
  manipulated_parnames <- character(length(context$manipulations))

  manipulated_parnames <- vapply(
    context$manipulations,
    function(x) as.character(x[[2]]),
    character(1)
  )
  names(context$manipulations) <- manipulated_parnames

  estimated_parnames <- setdiff(parnames, c(manipulated_parnames, fixed_parnames))

  # No overlaps
  if (length(intersect(estimated_parnames, manipulated_parnames)) > 0) {
    stop(sprintf(
      "Invalid parameter overlap: estimated and manipulated parameters overlap: %s",
      paste(intersect(estimated_parnames, manipulated_parnames), collapse = ", ")
    ))
  }

  if (length(intersect(estimated_parnames, fixed_parnames)) > 0) {
    stop(sprintf(
      "Invalid parameter overlap: estimated and fixed parnames overlap: %s",
      paste(intersect(estimated_parnames, fixed_parnames), collapse = ", ")
    ))
  }

  if (length(intersect(fixed_parnames, manipulated_parnames)) > 0) {
    stop(sprintf(
      "Invalid parameter overlap: fixed and manipulated parameters overlap: %s",
      paste(intersect(fixed_parnames, manipulated_parnames), collapse = ", ")
    ))
  }

  # No coverage mismatches
  all_assigned <- c(fixed_parnames, manipulated_parnames, estimated_parnames)
  if (!setequal(parnames, all_assigned)) {
    missing <- setdiff(parnames, all_assigned)
    extra <- setdiff(all_assigned, parnames)
    error_msg <- "Parameter coverage mismatch:"
    if (length(missing) > 0) error_msg <- paste0(error_msg, "\n- Missing: ", paste(missing, collapse = ", "))
    if (length(extra) > 0) error_msg <- paste0(error_msg, "\n- Unexpected: ", paste(extra, collapse = ", "))
    stop(error_msg)
  }

  # thetas_parnames
  base_names <- if (context$sym_thetas) "theta" else c("thetaLower", "thetaUpper")
  if (context$n_ratings <= 2) {
    thetas_parnames <- base_names
  } else {
    first_thetas <- paste0(base_names, "1")
    dtheta_indices <- 2:(context$n_ratings - 1)
    dthetas <- if (context$sym_thetas) {
      paste0("dtheta", dtheta_indices)
    } else {
      paste0(
        "dtheta",
        rep(c("Lower", "Upper"), times = context$n_ratings - 2),
        rep(dtheta_indices, times = 2)
      )
    }
    thetas_parnames <- c(first_thetas, dthetas)
  }

  context$beta_names <- c(estimated_parnames, thetas_parnames)
  context$estimated_params <- setNames(
    lapply(estimated_parnames, function(p) which(context$beta_names == p)),
    estimated_parnames
  )

  # Fill in default optimizer options if missing
  missing_opts <- !(names(DEFAULT_OPTS) %in% names(context$opts))
  context$opts <- c(context$opts, DEFAULT_OPTS[missing_opts])

  # simult_conf and restr_tau
  context$simult_conf <- FALSE
  if (context$restr_tau == "simult_conf") {
    context$simult_conf <- TRUE
    context$restr_tau <- 1
  }

  context
}
