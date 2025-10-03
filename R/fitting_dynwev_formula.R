#' @keywords internal
fitting_dynwev_formula <- function(context) {
  #### grid search setup ####
  if (context$grid_search) {
    n_initials <- 100
    inits <- matrix(
      rnorm(n_initials * length(context$beta_names)),
      nrow = n_initials,
      dimnames = list(NULL, context$beta_names)
    )
  } else {
    single_start <- rnorm(length(context$beta_names), mean = 0, sd = 0.1)
    inits <- matrix(
      single_start,
      nrow = 1,
      dimnames = list(NULL, context$beta_names)
    )
  }

  # initial setup logging
  if (context$logging) log_model_info(context)

  if ("st0" %in% context$beta_names) {
    if (context$logging) log_progress("Rescaling st0 initials to improve integration time")
    inits[, "st0"] <- inits[, "st0"] / 2 - 1.6
  }

  #### parallel cluster setup ####
  if (context$parallel) {
    cl <- makeCluster(context$n_cores, type = "SOCK")
    on.exit(try(stopCluster(cl), silent = TRUE))
    clusterExport(cl, varlist = c("context"), envir = environment())
    if (context$logging) log_progress(sprintf("Initialized parallel cluster with %d cores", context$n_cores))
  }

  #### grid search ####
  log_likelihood <- NULL
  if (context$grid_search) {
    if (context$logging) {
      log_progress(sprintf(
        "%d parameter sets to check for %d rows of data",
        nrow(inits), nrow(context$dependent_vars)
      ))
      log_progress("Starting grid search...")
    }

    inits_rows <- lapply(seq_len(nrow(inits)), function(i) inits[i, ])
    log_likelihood <- if (context$parallel) {
      parLapply(cl, inits_rows, function(row) grid_search_worker(context, row))
    } else {
      lapply(inits_rows, function(row) grid_search_worker(context, row))
    }

    log_likelihood <- vapply(log_likelihood, identity, numeric(1))
    inits <- inits[order(log_likelihood), ]

    if (context$logging) {
      log_progress(sprintf(
        "Grid search complete - Best likelihood: %s",
        format_numeric(min(log_likelihood))
      ), "success")

      save_optimization_state(context,
        list(log_likelihood = log_likelihood, inits = inits),
        type = "grid_search"
      )
    }
  }

  #### optimization ####
  if (context$logging) {
    log_progress(sprintf(
      "Starting optimization with %d attempts",
      context$opts$n_attempts
    ))
  }

  starts <- inits[seq_len(context$opts$n_attempts), , drop = FALSE]
  starts_rows <- lapply(seq_len(nrow(starts)), function(i) starts[i, ])

  optim_outs <- if (context$parallel && context$opts$n_attempts > 1) {
    clusterExport(cl, varlist = c("optim_node"), envir = environment())
    parLapply(
      cl,
      starts_rows,
      function(start_params) optim_node(context, start_params)
    )
  } else {
    lapply(
      starts_rows,
      function(start_params) optim_node(context, start_params)
    )
  }

  values <- vapply(optim_outs, function(x) x$value, numeric(1))
  best_idx <- which.min(values)
  fit <- if (length(best_idx) == 0) NULL else optim_outs[[best_idx]]

  if (context$logging) {
    save_optimization_state(context,
      list(
        optim_outs = optim_outs,
        best_idx = best_idx,
        values = values
      ),
      type = "optimization"
    )
  }

  #### wrap up results ####
  res <- list()
  if (!is.null(fit) && !is.na(fit$value)) {
    if (!is.null(context$used_ratings)) {
      fit$par <- fill_thresholds(fit$par, context$used_ratings, context$initial_n_ratings, context$sym_thetas)
    }

    res <- build_result_list(context, fit)

    if (context$logging) {
      log_progress("Optimization complete", "success")
      log_progress(sprintf("Final negLogLik: %s", format_numeric(res$negLogLik)))
      log_progress(sprintf("BIC: %s", format_numeric(res$BIC)))
      save_optimization_state(context,
        list(fit = fit, results = res),
        type = "final"
      )
    }
  } else {
    if (context$logging) {
      log_progress("No valid fit obtained - all attempts returned NA", "warning")
    } else {
      warning("No valid fit could be obtained, all optimization attempts returned NA")
    }
  }

  res
}

#' @keywords internal
optim_node <- function(context, start_params) {
  if (context$logging) {
    log_progress(sprintf(
      "Starting optimization node with %d restarts - Initial params: %s",
      context$opts$n_restarts,
      paste(names(start_params), format_numeric(start_params),
        sep = "=", collapse = ", "
      )
    ))
  }

  node_fit <- NULL
  best_value <- Inf

  for (j in seq_len(context$opts$n_restarts)) {
    # Jitter parameters
    jittered_params <- start_params +
      rnorm(length(start_params), sd = pmax(0.001, abs(start_params / 20)))

    if (context$logging) {
      log_progress(sprintf(
        "Restart %d/%d - Jittered params: %s",
        j, context$opts$n_restarts,
        paste(
          names(jittered_params),
          format_numeric(jittered_params),
          sep = "=", collapse = ", "
        )
      ))
    }

    current_fit <- tryCatch({
      nlopt_optimizer(context, jittered_params)
    }, error = function(e) {
      if (context$logging) {
        log_progress(sprintf(
          "Optimization failed on restart %d: %s",
          j, e$message
        ), "error")
      }
      NULL
    })

    if (!is.null(current_fit) && !is.na(current_fit$value) && current_fit$value < best_value) {
      node_fit <- current_fit
      best_value <- current_fit$value

      if (context$logging) {
        log_progress(sprintf(
          "New best fit on restart %d - negLogLik: %s - Params: %s",
          j, format_numeric(best_value),
          paste(names(current_fit$par), format_numeric(current_fit$par),
            sep = "=", collapse = ", "
          )
        ), "success")
      }
    }

    if (!is.null(node_fit)) {
      start_params <- node_fit$par
    }
  }

  if (is.null(node_fit)) {
    if (context$logging) {
      log_progress("Node failed to find valid fit across all restarts", "warning")
    }
    list(
      value = NA_real_,
      par = setNames(rep(NA_real_, length(start_params)), names(start_params))
    )
  } else {
    if (context$logging) {
      log_progress(sprintf(
        "Node complete - Final negLogLik: %s",
        format_numeric(node_fit$value)
      ))
    }
    node_fit
  }
}

#' @keywords internal
build_result_list <- function(context, fit) {
  k <- length(fit$par)
  N <- nrow(context$dependent_vars)
  AIC <- 2 * fit$value + 2 * k

  list(
    k = k,
    N = N,
    beta = fit$par,
    fixed = paste(
      c("sym_thetas", names(context$fixed)),
      c(context$sym_thetas, unlist(context$fixed)),
      sep = "=", collapse = ", "
    ),
    negLogLik = fit$value,
    BIC = 2 * fit$value + k * log(N),
    AIC = AIC,
    AICc = AIC + (2 * k * (k + 1)) / (N - k - 1)
  )
}
