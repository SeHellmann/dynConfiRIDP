#' @keywords internal
#' @importFrom stats
#'   rnorm
#'   setNames
#' @importFrom parallel
#'   makeCluster
#'   stopCluster
#'   clusterExport
#'   parLapply
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

  # Rescale st0 initials to lower values, because integration would otherwise take a lot of time
  if ("st0" %in% context$beta_names) inits[, "st0"] <- inits[, "st0"] / 2 - 1.6

  #### parallel cluster setup ####
  if (context$parallel) {
    cl <- makeCluster(context$n_cores, type = "SOCK")
    on.exit(try(stopCluster(cl), silent = TRUE))
    clusterExport(cl, varlist = c("context"), envir = environment())
  }

  #### grid search ####
  log_likelihood <- NULL
  if (context$grid_search) {
    if (context$logging) {
      logger::log_info(sprintf(
        "%d parameter sets to check for %d rows of data",
        nrow(inits),
        nrow(context$dependent_vars)
      ))
      logger::log_info("Searching initial values...")
      start_time <- Sys.time()
    }
    inits_rows <- lapply(seq_len(nrow(inits)), function(i) inits[i, ])

    log_likelihood <- if (context$parallel) {
      parallel::parLapply(cl, inits_rows, function(row) grid_search_worker(context, row))
    } else {
      lapply(inits_rows, function(row) grid_search_worker(context, row))
    }

    if (context$logging) {
      logger::log_success(sprintf(
        "Intial grid search took %.2f",
        difftime(Sys.time(), start_time, units = "mins")
      ))
      save(log_likelihood, inits, context$dependent_vars, file = context$logfile)
    }

    log_likelihood <- vapply(log_likelihood, identity, numeric(1))
    inits <- inits[order(log_likelihood), ]
  }

  #### optimization ####
  if (context$logging) logger::log_info("Start fitting...")

  starts <- inits[seq_len(context$opts$n_attempts), , drop = FALSE]
  starts_rows <- lapply(seq_len(nrow(starts)), function(i) starts[i, ])

  optim_outs <- if (context$parallel && context$opts$n_attempts > 1) {
    clusterExport(cl, varlist = c("optim_node"), envir = environment())
    parallel::parLapply(
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
    save(log_likelihood, context$dependent_vars, inits, optim_outs, file = context$logfile)
  }

  #### wrap up results ####
  res <- list()
  if (!is.null(fit) && !is.na(fit$value)) {
    if (!is.null(context$used_ratings)) {
      fit$par <- fill_thresholds(fit$par, context$used_ratings, context$initial_n_ratings, context$sym_thetas)
    }

    res$k <- length(fit$par)
    res$N <- nrow(context$dependent_vars)
    res$beta <- fit$par
    res$fixed <- paste(
      c("sym_thetas", names(context$fixed)),
      c(context$sym_thetas, unlist(context$fixed)),
      sep = "=", collapse = ", "
    )
    res$negLogLik <- fit$value
    res$BIC <- 2 * res$negLogLik + res$k * log(res$N)
    res$AIC <- 2 * res$negLogLik + 2 * res$k
    res$AICc <- res$AIC + (2 * res$k * (res$k + 1)) / (res$N - res$k - 1)

    if (context$logging) {
      logger::log_success("Done fitting and autosaved results")
      save(log_likelihood, context$dependent_vars, context$model_matrix, fit, inits, res, file = context$logfile)
    }
  } else {
    if (context$logging) {
      logger::log_warn("No valid fit could be obtained, all optimization attempts returned NA")
    } else {
      warning("No valid fit could be obtained, all optimization attempts returned NA")
    }
  }

  res
}

#' @keywords internal
optim_node <- function(context, start_params) {
  node_fit <- NULL
  for (j in seq_len(context$n_restarts)) {
    # jitter start
    start_params <- start_params + rnorm(length(start_params), sd = pmax(0.001, abs(start_params / 20)))
    # STUB for nlopt rcpp export
    m <- tryCatch({
      nlopt_optimizer(context, start_params)
    }, error = function(e) {
      if (context$logging) logger::log_error(paste("Optimization failed:", e$message))
      NULL
    })

    if (!is.null(m) && (is.null(node_fit) || m$value < node_fit$value)) {
      node_fit <- m
    }

    if (!is.null(node_fit)) start_params <- node_fit$par
  }

  if (is.null(node_fit)) {
    list(
      value = NA_real_,
      par = setNames(rep(NA_real_, length(start_params)), names(start_params))
    )
  } else {
    node_fit
  }
}