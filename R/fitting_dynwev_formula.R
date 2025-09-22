#' @keywords internal
#' @importFrom stats
#'   rnorm
#' @importFrom parallel
#'   makeCluster
#'   stopCluster
#'   clusterExport
#'   parApply
fitting_dynwev_formula <- function(context) {
  beta_names <- c(context$fit_beta_parnames, context$const_parnames, context$thetas_parnames)
  n_initials <- 100
  inits <- matrix(rnorm(n_initials * length(beta_names)), nrow = n_initials, dimnames = list(NULL, beta_names))
  # Rescale st0 initials to lower values, because integration would otherwise take a lot of time
  if ("st0" %in% beta_names) inits[, "st0"] <- inits[, "st0"] / 2 - 1.6
  if (!context$grid_search) inits <- colMeans(inits)

  #### setup cluster for parallelization ####
  cl <- NULL
  if (context$parallel) {
    cl <- makeCluster(context$n_cores, type = "PSOCK")
    on.exit(try(stopCluster(cl), silent = TRUE))
  }

  #### grid search ####
  log_likelihood <- NULL
  if (context$parallel) {
    if (context$logging) {
      logger::log_info(sprintf(
        "%d parameter sets to check for %d rows of data",
        nrow(inits),
        nrow(context$dependent_vars)
      ))
      logger::log_info("Searching initial values...")
      start_time <- Sys.time()
    }

    # STUBS for neglikelihood_formula()
    log_likelihood <- if (context$parallel) {
      parApply(cl, inits, MARGIN = 1, function() {})
    } else {
      apply(inits, MARGIN = 1, function() {})
    }

    if (context$logging) {
      logger::log_success(sprintf(
        "Intial grid search took %.2f",
        difftime(Sys.time(), start_time, units = "mins")
      ))
      save(log_likelihood, inits, context$dependent_vars, file = context$logfile)
    }

    log_likelihood <- as.numeric(log_likelihood)
    inits <- inits[order(log_likelihood), ]
  }

  #### optimization ####
  if (context$logging) logger::log_info("Start fitting...")

  fit <- NULL
  if (!context$parallel || (context$opts$n_attempts == 1)) {
    for (i in seq_len(context$opts$n_attempts)) {
      start <- inits[i, ]
      for (j in seq_len(context$opts$n_restarts)) {
        # jitter start
        start <- start + rnorm(length(start), sd = pmax(0.001, abs(start / 20)))
        # STUB for nlopt rcpp export
        m <- tryCatch(function() {}, error = function(e) NULL)

        if (!is.null(m) && (is.null(fit) || m$value < fit$value)) {
          fit <- m
          if (context$logging) {
            logger::log_info(sprintf("New best fit at attempt %d, restart %d", i, j))
            save(log_likelihood, inits, context$dependent_vars, fit, file = context$logfile)
          }
        }
        if (!is.null(fit)) start <- fit$par
      }
    }
  } else {
    starts <- inits[seq_len(context$opts$n_attempts), , drop = FALSE]

    optim_node <- function(start_params) {
      node_fit <- NULL
      for (j in seq_len(context$opts$n_restarts)) {
        # jitter start
        start_params <- start_params + rnorm(length(start_params), sd = pmax(0.001, abs(start_params / 20)))
        # STUB for nlopt rcpp export
        m <- tryCatch(function() {}, error = function(e) NULL)

        if (!is.null(m) && (is.null(node_fit) || m$value < node_fit$value)) {
          node_fit <- m
        }
        if (!is.null(node_fit)) start_params <- node_fit$par
      }
      if (is.null(node_fit)) c(NA, rep(NA, length(start_params))) else c(node_fit$value, node_fit$par)
    }

    optim_outs <- parApply(cl, starts, MARGIN = 1, optim_node)
    stopCluster(cl)

    save(log_likelihood, context$dependent_vars, inits, optim_outs, file = context$logfile)

    optim_outs <- t(optim_outs)
    best_res <- optim_outs[which.min(optim_outs[, 1]), ]
    fit <- list(par = best_res[-1], value = best_res[1])
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
  }

  res
}
