#' @keywords internal
setup_parallel <- function(context) {
  if (!context$parallel) return(context)

  n_cores <- if (is.null(context$n_cores)) availableCores() - 1 else context$n_cores
  if (future::supportsMulticore()) {
    plan_type <- "multicore"
    future::plan(multicore, workers = n_cores)
  } else {
    plan_type <- "multisession"
    future::plan(multisession, workers = n_cores)
  }

  if (context$logging) logger::log_info(sprintf("Initializing %s plan with %d workers", plan_type, n_cores))

  context
}