#' @keywords internal
setup_parallel <- function(context) {
  if (!context$parallel) return(context)

  n_available <- availableCores()
  n_requested <- 1

  n_cores <- if (is.null(context$n_cores)) n_available - 1 else context$n_cores

  if (length(n_cores) > 1) {
    #### nested parallel setup ####
    n_outer <- n_cores[1]
    n_inner <- n_cores[2]
    n_requested <- n_outer * n_inner

    plan_strategy <- list(
      tweak(multisession, workers = n_outer),
      tweak(multisession, workers = n_inner)
    )
    plan(plan_strategy)

    zombie_msg <- paste(
      "Pressing Ctrl-C will NOT terminate the parallel workers.",
      "They will become 'zombie' processes and must be manually killed."
    )

    if (context$logging) {
      log_info(sprintf(
        "Initializing NESTED parallel plan: %d outer workers, %d inner workers",
        n_outer, n_inner
      ))
      log_warn(zombie_msg)
    } else {
      warning(zombie_msg)
    }
  } else {
    #### flat parallel setup ####
    n_single <- n_cores[1]
    n_requested <- n_single

    if (supportsMulticore()) {
      plan_type <- "multicore"
      plan(multicore, workers = n_single)
    } else {
      plan_type <- "multisession"
      plan(multisession, workers = n_single)
    }

    if (context$logging) {
      log_info(sprintf(
        "Initializing FLAT %s plan with %d workers",
        plan_type, n_single
      ))
    }
  }

  if (n_requested > n_available) {
    warn_msg <- sprintf(
      "A total of %d parallel workers were requested, but only %d cores are available.",
      "This may lead to over-subscription and poor performance.",
      n_requested, n_available
    )

    if (context$logging) {
      log_warn(warn_msg)
    } else {
      warning(warn_msg)
    }
  }

  context
}