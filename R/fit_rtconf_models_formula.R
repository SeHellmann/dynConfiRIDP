#' Fit Multiple Sequential Sampling Confidence Models
#'
#' A wrapper function to fit multiple models for multiple participants,
#' typically in parallel.
#'
#' @details
#' This function fits sequential sampling models to data from multiple
#' participants. It splits the data by participant ID and calls
#' `fit_rtconf_formula` for each participant-model combination.
#'
#' **Data Requirements:**
#' The input `data` must be a `data.frame` with columns for:
#' * **Subject ID**: One of `sbj`, `participant`, or `subject`.
#' * `rt`: Response time (a positive numeric vector).
#' * `rating`: Discrete confidence rating (an integer vector >= 0, with at
#'     least 2 unique levels).
#' * **Decision Outcome**: One of the following combinations:
#'     * `response` and `stimulus`: Columns indicating the participant's
#'         response and the correct stimulus. Each must have exactly
#'         two unique values.
#'     * `correct`: A column with `0` (error) and `1` (correct).
#'         If `stimulus` is also provided, `response` is inferred. If
#'         `stimulus` is *not* provided, bias parameters (like `z`)
#'         may not be identifiable.
#'
#' **Parallel Execution:**
#' If `parallel = TRUE`, fitting is parallelized using the `future` framework.
#' * **Flat Plan (`n_cores = 8`):** If `n_cores` is a single integer (>= 1),
#'     fitting is parallelized *across* the model-subject combinations.
#' * **Nested Plan (`n_cores = c(4, 5)`):** If `n_cores` is a vector of
#'     two integers (>= 1), a nested plan is created using
#'     `future::tweak(multisession, ...)` for both levels.
#'     * The **outer level** (`n_outer = 4`) runs model-subject
#'         combinations in parallel.
#'     * The **inner level** (`n_inner = 5`) parallelizes the grid_search and the
#'         `opts$n_attempts` (different starting values) *within*
#'         each model-subject fit. This value is ideally set to
#'         match your `opts$n_attempts`.
#'     * The `opts$n_restarts` for each attempt will still run sequentially.
#'
#' **Manipulations:**
#' The `manipulations` argument is a list of formulas that defines
#' experimental manipulations It can be provided in two ways:
#' 1.  A single list of formulas (e.g., `list(v ~ condition)`). This *same*
#'     list will be applied to *all* models.
#' 2.  A list of lists, where each inner list corresponds to a model in the
#'     `models` argument (e.g., `list(list(v ~ cond1), list(v ~ cond2))`).
#'     This must be the same length as `models`.
#'
#' @param data A `data.frame` with trial-by-trial data for one or more
#'   participants. See the "Data Requirements" section in `@details`.
#' @param models A character vector of model names to fit.
#'   Default: `c("dynaViTE", "dynWEV")`.
#'   Supported models are: "2DSD", "dynWEV", and "dynaViTE".
#' @param optim_method The optimization algorithm. Supports "Nelder-Mead"
#'   (default) and "bobyqa".
#' @param fixed A named `list` of model parameters to hold constant
#'   *during fitting*, where each value is a single numeric scalar
#'   (e.g., `list(s = 1, z = 0.5)`).
#' @param manipulations A `list` of formulas specifying how model parameters
#'   depend on predictors. See the "Manipulations" section in `@details`.
#' @param n_ratings The total number of distinct confidence ratings. Must be
#'   an integer >= 2. If `NULL` (default), this is inferred from the `rating`
#'   column.
#' @param restr_tau A constraint on the `tau` parameter. Must be a positive
#'   numeric value (e.g., `Inf`) or the string "simult_conf".
#'   Default: `Inf`.
#' @param sym_thetas Logical. If `TRUE`, confidence thresholds are assumed
#'   to be symmetric for both responses. Default: `FALSE`.
#' @param precision Numeric. An integer (e.g., 1, 2, 3...) that sets the
#'   accuracy level for the internal density calculations. Higher values
#'   increase accuracy and computation time, while lower values are faster
#'   but less precise. Default: 2.
#' @param opts A `list` controlling the optimization process:
#'   * `n_attempts` (integer >= 1): Number of optimization attempts.
#'       Default: 5.
#'   * `n_restarts` (integer >= 1): Number of restarts for each attempt.
#'       Default: 5.
#'   * `maxfun` (integer >= 1): Max function evaluations for the `nlopt` optimization.
#'       Default: 8000.
#'   * `reltol` (numeric > 0): Relative tolerance for the `nlopt` optimization.
#'       Default: 1e-6.
#' @param grid_search Logical. If `TRUE` (default), performs an initial
#'   grid search for each fit.
#' @param logging Logical. If `TRUE`, creates a `.logs` directory and
#'   saves detailed logs and intermediate results for each fit.
#' @param parallel Logical. If `TRUE`, parallelizes fitting across
#'   models and subjects.
#' @param n_cores Integer or vector. `NULL` (default) uses
#'   `future::availableCores() - 1`. See the "Parallel Execution" section
#'   in `@details` for flat vs. nested plans.
#'
#' @return A `list` where each element is the result from
#'   `fit_rtconf_formula` for a single participant-model combination.
#'
#' @examples
#' # 1. Generate data from two artificial participants (using SATdata)
#' \dontrun{
#' data <- SATdata
#' # Add a dummy participant column for the example
#' data$participant <- rep(1:2, each = nrow(data)/2)
#' data$condition <- factor(data$condition)
#'
#' # 2. Define manipulations
#' manipulations <- list(v ~ condition, a ~ SAT)
#'
#' # 3. Fit models
#' # This will fit "dynaViTE" and "2DSD" for both participant 1 and 2
#' # using 4 parallel workers.
#'
#' fixed_params <- list(z = 0.5, svis = 1, d = 0)
#'
#' fit_list <- fit_rtconf_models_formula(
#'   data,
#'   models = c("dynaViTE", "2DSD"),
#'   manipulations = manipulations,
#'   fixed = fixed_params,
#'   n_ratings = 6,
#'   sym_thetas = TRUE,
#'   logging = FALSE,
#'   parallel = TRUE,
#'   n_cores = 4
#' )
#'
#' # View the results for the first fit
#' print(fit_list[[1]])
#' }
#'
#' @export
#' @md
#' @rdname fit_rtconf_models_formula
#' @author Sebastian Hellmann, Manuel Methasani
#' @references
#' Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling
#' of choice, confidence, and response time in visual perception.
#' *Psychological Review*. Advance online publication.
#' <https://doi.org/10.1037/rev0000411>
fit_rtconf_models_formula <- function(
  data,
  models = c("dynaViTE", "dynWEV"),
  optim_method = "Nelder-Mead",
  fixed = list("s" = 1),
  manipulations = list(),
  n_ratings = NULL,
  restr_tau = Inf,
  sym_thetas = FALSE,
  precision = 2,
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
    setup_main_logging() |>
    setup_parallel()

  jobs_setup_res <- setup_jobs(base_context, models, manipulations)
  base_context <- jobs_setup_res$context
  jobs_list <- jobs_setup_res$jobs_list

  if (base_context$parallel) on.exit(plan(sequential), add = TRUE)

  future_lapply(
    jobs_list,
    function(job) fit_rtconf_formula_worker(base_context, job),
    future.seed = TRUE,
    future.scheduling = FALSE
  )
}

#' @keywords internal
fit_rtconf_formula_worker <- function(base_context, job) {
  job_context <- base_context

  job_context$data <- subset(base_context$data, base_context$data[["sbj"]] == job$sbj)
  job_context$model <- job$model
  job_context$manipulations <- job$manipulations

  job_context <- setup_subject_logging(job_context)

  res <- fit_rtconf_formula_dispatcher(job_context)

  res$model <- job$model
  res$sbj <- job$sbj

  res
}
