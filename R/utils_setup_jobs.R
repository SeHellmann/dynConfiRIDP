#' @keywords internal
setup_jobs <- function(context, models, manipulations) {
  # determine number of jobs, i.e. model-participant-combinations
  potential_sbj_cols <- c("sbj", "participant", "subject")
  sbjcol <- potential_sbj_cols[potential_sbj_cols %in% names(context$data)][1]
  if (is.na(sbjcol)) {
    context$data$sbj <- 999
  } else {
    if (sbjcol != "sbj") {
      context$data$sbj <- context$data[[sbjcol]]
      context$data[[sbjcol]] <- NULL
    }
  }
  subjects <- unique(context$data$sbj)

  # prepare the job list for combinations of models and subjects
  jobs <- expand.grid(model = seq_along(models), sbj = subjects)
  jobs_list <- vector("list", nrow(jobs))

  for (i in seq_len(nrow(jobs))) {
    model_idx <- jobs$model[i]

    # handle the two ways manipulations can be provided:
    # 1. single list of formulas (applied to all models)
    if (length(manipulations) == 0 || inherits(manipulations[[1]], "formula")) {
      curr_manipulations <- manipulations
    } else {
      # 2. list-of-lists (one for each model)
      curr_manipulations <- manipulations[[model_idx]]
    }
    jobs_list[[i]] <- list(
      model = models[model_idx],
      sbj = jobs$sbj[i],
      manipulations = curr_manipulations
    )
  }

  list(
    context = context,
    jobs_list = jobs_list
  )
}