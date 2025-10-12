#' @keywords internal
setup_parallel <- function(context) {
  if (!context$parallel) return(context)

  n_cores <- if (is.null(context$n_cores)) availableCores() - 1 else context$n_cores
  if (supportsMulticore()) {
    plan_type <- "multicore"
    plan(multicore, workers = n_cores)
  } else {
    plan_type <- "multisession"
    plan(multisession, workers = n_cores)
  }

  if (context$logging) log_info(sprintf("Initializing %s plan with %d workers", plan_type, n_cores))

  context
}