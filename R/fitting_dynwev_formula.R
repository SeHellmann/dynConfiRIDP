#' @keywords internal
fitting_dynwev_formula <- function(context) {
  optimization_context <- get_dynwev_optimization_context(context)

  #### grid search setup ####
  if (context$grid_search) {
    n_initials <- 100
    inits <- matrix(
      rnorm(n_initials * length(context$beta_names)),
      nrow = n_initials,
      dimnames = list(NULL, context$beta_names)
    )
  } else {
    inits <- matrix(
      rnorm(context$opts$n_attempts * length(context$beta_names), mean = 0, sd = 0.1),
      nrow = context$opts$n_attempts,
      dimnames = list(NULL, context$beta_names)
    )
  }

  # initial setup logging
  if (context$logging) {
    log_info(sprintf("Model: %s", context$model))
    log_info(sprintf("Data rows: %d", nrow(context$data)))
  }

  if ("st0" %in% context$beta_names) {
    if (context$logging) log_info("Rescaling st0 initials to improve integration time")
    inits[, "st0"] <- inits[, "st0"] / 2 - 1.6
  }

  #### grid search ####
  log_likelihood <- NULL
  if (context$grid_search) {
    if (context$logging) {
      log_info(sprintf(
        "%d parameter sets to check for %d rows of data",
        n_initials, nrow(context$dependent_vars)
      ))
      log_info("Starting grid search...")
    }

    min_batches <- 10
    n_workers <- future::nbrOfWorkers()
    n_batches <- max(n_workers, min_batches)

    batch_row_indices <- split(
      seq_len(n_initials),
      cut(seq_len(n_initials), breaks = n_batches, labels = FALSE)
    )

    inits_batches <- lapply(
      batch_row_indices,
      function(indices) inits[indices, , drop = FALSE]
    )

    log_likelihood_list <- future_lapply(
      inits_batches,
      function(batch) grid_search_worker(optimization_context, batch),
      future.seed = TRUE
    )

    log_likelihood <- unlist(log_likelihood_list, use.names = FALSE)
    inits <- inits[order(log_likelihood), ]

    if (context$logging) {
      log_success(sprintf(
        "Grid search complete - Best negLogLik: %.4f",
        min(log_likelihood)
      ))

      save_optimization_state(context,
        list(
          context = context,
          log_likelihood = log_likelihood,
          inits = inits
        ),
        type = "grid_search"
      )
    }
  }

  #### optimization ####
  if (context$logging) {
    log_info(sprintf(
      "Starting optimization with %d attempts",
      context$opts$n_attempts
    ))
  }

  starts <- inits[seq_len(context$opts$n_attempts), , drop = FALSE]
  starts_rows <- lapply(seq_len(nrow(starts)), function(i) starts[i, ])
  names(starts_rows) <- paste0("attempt_", seq_along(starts_rows))

  optim_outs <- future_lapply(
    starts_rows,
    function(start_params) optimization_node(optimization_context, start_params),
    future.seed = TRUE
  )

  node_values <- vapply(optim_outs, function(node_result) {
    if (is.na(node_result$best_fit_idx)) return(Inf)

    best_fit <- node_result$all_fits[[node_result$best_fit_idx]]
    best_fit$value
  }, numeric(1))

  best_node_idx <- which.min(node_values)

  if (context$logging) {
    save_optimization_state(context,
      list(
        optim_outs = optim_outs,
        best_node_idx = best_node_idx,
        node_values = node_values
      ),
      type = "optimization"
    )
  }

  best_node_result <- optim_outs[[best_node_idx]]
  fit <- best_node_result$all_fits[[best_node_result$best_fit_idx]]

  #### wrap up results ####
  res <- list()
  if (!is.na(fit$value) && !is.null(fit$params)) {
    if (!is.null(context$used_ratings)) {
      fit$params <- fill_thresholds(fit$params, context$used_ratings, context$initial_n_ratings, context$sym_thetas)
    }

    res <- build_result_list(context, fit)

    if (context$logging) {
      log_success("Optimization complete")
      log_info(sprintf("Final negLogLik: %.4f", res$negLogLik))
      log_info(sprintf("BIC: %.4f", res$BIC))
      save_optimization_state(context,
        list(fit = fit, results = res),
        type = "final"
      )
    }
  } else {
    if (context$logging) {
      log_warn("No valid fit could be obtained, all optimization attempts returned NA")
    } else {
      warning("No valid fit could be obtained, all optimization attempts returned NA")
    }
  }

  res
}

#' @keywords internal
optimization_node <- function(optimization_context, start_params) {
  if (optimization_context$logging) {
    log_info(sprintf(
      "Starting optimization node with %d restarts - Initial params: %s",
      optimization_context$opts$n_restarts,
      paste(names(start_params), sprintf("%.4f", start_params),
        sep = "=", collapse = ", "
      )
    ))
  }

  all_fits <- vector("list", optimization_context$opts$n_restarts)
  best_fit_idx <- NA_integer_
  current_params <- start_params

  for (i in seq_len(optimization_context$opts$n_restarts)) {
    # Jitter parameters
    jittered_params <- current_params +
      rnorm(length(current_params), sd = pmax(0.001, abs(current_params / 20)))

    if (optimization_context$logging) {
      log_info(sprintf(
        "Restart %d/%d - Jittered params: %s",
        i, optimization_context$opts$n_restarts,
        paste(
          names(jittered_params),
          sprintf("%.4f", jittered_params),
          sep = "=", collapse = ", "
        )
      ))
    }

    current_fit <- nlopt_optimizer(optimization_context, jittered_params)
    all_fits[[i]] <- current_fit

    if (!is.na(current_fit$value) && !is.null(current_fit$params)) {
      if (is.na(best_fit_idx) || current_fit$value < all_fits[[best_fit_idx]]$value) {
        current_params <- current_fit$params
        best_fit_idx <- i

        if (optimization_context$logging) {
          log_success(sprintf(
            "New best fit on restart %d - negLogLik: %s - Params: %s",
            i, sprintf("%.4f", current_fit$value),
            paste(names(current_fit$params), sprintf("%.4f", current_fit$params),
              sep = "=", collapse = ", "
            )
          ))
        }
      }
    }
  }

  if (optimization_context$logging) {
    if (is.na(best_fit_idx)) {
      log_warn("Node failed to find valid fit across all restarts")
    } else {
      log_info(sprintf(
        "Node complete - Best negLogLik: %.4f",
        all_fits[[best_fit_idx]]$value
      ))
    }
  }

  list(
    all_fits = all_fits,
    best_fit_idx = best_fit_idx
  )
}

#' @keywords internal
build_result_list <- function(context, fit) {
  k <- length(fit$params)
  N <- nrow(context$dependent_vars)
  AIC <- 2 * fit$value + 2 * k

  list(
    k = k,
    N = N,
    beta = fit$params,
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
