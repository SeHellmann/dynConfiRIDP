#' @export
fit_rtconf_models_formula <- function(
  data,
  models = c("dynaViTE", "PCRMt", "IRMt"),
  optim_method = "Nelder-Mead",
  fixed = list("s" = 1),
  manipulations = list(),
  n_ratings = NULL,
  restr_tau = Inf,
  sym_thetas = FALSE,
  precision = 1e-5,
  opts = list(),
  grid_search = TRUE,
  logging = FALSE,
  parallel = FALSE,
  n_cores = NULL
) {
  validate_rtconf_models_args(
    models,
    optim_method,
    manipulations,
    parallel,
    n_cores
  )

  ### Determine number of jobs, i.e. model-participant-combinations
  potential_sbj_cols <- c("subject", "participant", "sbj")
  sbjcol <- potential_sbj_cols[match(TRUE, potential_sbj_cols %in% names(data))]
  if (length(sbjcol) == 0) {
    data$sbj <- 999
    sbjcol <- "sbj"
  } else {
    if (sbjcol != "sbj") {
      data$sbj <- data[[sbjcol]]
      data[[sbjcol]] <- NULL
    }
  }
  subjects <- unique(data$sbj)

  # Prepare the job list for combinations of models & subjects
  jobs <- expand.grid(model = seq_along(models), sbj = subjects)
  jobs_list <- vector("list", nrow(jobs))

  for (i in seq_len(nrow(jobs))) {
    model_idx <- jobs$model[i]
    if (length(manipulations) == 0 || inherits(manipulations[[1]], "formula")) {
      curr_manipulations <- manipulations
    } else {
      curr_manipulations <- manipulations[[model_idx]]
    }
    jobs_list[[i]] <- list(
      model = models[model_idx],
      sbj = jobs$sbj[i],
      manipulations = curr_manipulations
    )
  }

  call_fitfct <- function(job) {
    sbj_data <- subset(data, data[["sbj"]] == job$sbj)

    res <- fit_rtconf_formula(
      data = sbj_data,
      model = job$model,
      optim_method = optim_method,
      fixed = fixed,
      manipulations = job$manipulations,
      n_ratings = n_ratings,
      restr_tau = restr_tau,
      sym_thetas = sym_thetas,
      precision = precision,
      opts = opts,
      grid_search = grid_search,
      logging = logging
    )
    res$model <- job$model
    res$sbj <- job$sbj

    res
  }

  # TODO: integrate better with fit_rtconf_formula setup
  context <- list(
    logging = logging,
    parallel = parallel,
    n_cores = n_cores
  )

  setup_logging(context)
  setup_parallel(context)


  if (parallel) on.exit(plan(sequential), add = TRUE)

  res <- future_lapply(
    jobs_list,
    function(job) call_fitfct(job),
    future.seed = TRUE,
  )

  res
}
