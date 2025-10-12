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
  args <- as.list(environment())
  do.call(validate_rtconf_models_args, args)

  base_context <- args |>
    get_base_context() |>
    setup_logging() |>
    setup_parallel()

  jobs_setup_res <- setup_jobs(base_context, models, manipulations)
  base_context <- jobs_setup_res$context
  jobs_list <- jobs_setup_res$jobs_list

  if (base_context$parallel) on.exit(plan(sequential), add = TRUE)

  future_lapply(
    jobs_list,
    function(job) fit_rtconf_formula_worker(base_context, job),
    future.seed = TRUE,
  )
}

#' @keywords internal
fit_rtconf_formula_worker <- function(base_context, job) {
  job_context <- base_context

  job_context$data <- subset(base_context$data, base_context$data[["sbj"]] == job$sbj)
  job_context$model <- job$model
  job_context$manipulations <- job$manipulations

  res <- fit_rtconf_formula_dispatcher(job_context)

  res$model <- job$model
  res$sbj <- job$sbj

  res
}
