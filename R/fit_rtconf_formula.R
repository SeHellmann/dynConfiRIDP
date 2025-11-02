#' Fit a Sequential Sampling Confidence Model
#'
#' Fits the parameters of a single sequential sampling model for response time
#' and confidence, where experimental manipulations are specified using formulas.
#'
#' @details
#' This function provides a formula-based interface for fitting dynamic models of
#' decision confidence. It supports models like dynWEV, dynaViTE, and 2DSD.
#' The fitting process involves an initial grid search
#' (if `grid_search = TRUE`) followed by an optimization using `nlopt`.
#'
#' **Data Requirements:**
#' The input `data` must be a `data.frame` with columns for:
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
#' **Parameter Manipulations:**
#' The `manipulations` argument allows parameters to vary across conditions.
#' It takes a list of formulas, e.g., `list(v ~ condition, a ~ stimulus)`.
#' The fitting function uses `stats::model.matrix` to build the predictor
#' matrix.
#'
#' **Parallel Execution:**
#' If `parallel = TRUE`, fitting is parallelized using the `future` framework.
#' The parallelization happens over the `opts$n_attempts` (different starting
#' values). The `opts$n_restarts` for each attempt still run sequentially.
#'
#' @param data A `data.frame` with trial-by-trial data. See the
#'   "Data Requirements" section in `@details`.
#' @param model A character string specifying the model to fit.
#'   Default: "dynWEV".
#'   Supported models are: "2DSD", "dynWEV", and "dynaViTE".
#' @param optim_method The optimization algorithm to use.
#'   Default: "Nelder-Mead". Supports "Nelder-Mead" and "bobyqa".
#' @param fixed A named `list` of model parameters to hold constant
#'   *during fitting*, where each value is a single numeric scalar
#'   (e.g., `list(s = 1, z = 0.5)`). Default: `list("s" = 1)`.
#' @param manipulations A `list` of formulas specifying how model parameters
#'   depend on predictors in `data` (e.g., `list(v ~ condition)`). Default: `list()`.
#' @param n_ratings The total number of distinct confidence ratings. Must be
#'   an integer >= 2. If `NULL` (default), this is inferred from the `rating`
#'   column.
#' @param restr_tau A constraint on the `tau` parameter. Must be a positive
#'   numeric value (e.g., `Inf`) or the string "simult_conf".
#'   Default: `Inf`.
#' @param sym_thetas Logical. If `TRUE`, confidence thresholds are assumed
#'   to be symmetric. Default: `FALSE`.
#' @param precision Numeric. An integer (e.g., 1, 2, 3...) that sets the
#'   accuracy level for the internal density calculations. Higher values
#'   increase accuracy and computation time, while lower values are faster
#'   but less precise. Default: 2.
#' @param opts A `list` controlling the optimization process:
#'   * `n_attempts` (integer >= 1): Number of optimization attempts.
#'       Default: 5.
#'   * `n_restarts` (integer >= 1): Number of restarts for each attempt.
#'       Default: 5.
#'   * `maxfun` (integer >= 1): Max function evaluations for `nlopt`.
#'       Default: 8000.
#'   * `reltol` (numeric > 0): Relative tolerance for `nlopt`.
#'       Default: 1e-6.
#' @param grid_search Logical. If `TRUE` (default), performs an initial
#'   grid search to find good starting parameters.
#' @param logging Logical. If `TRUE`, creates a `.logs` directory and
#'   saves detailed logs and intermediate results.
#' @param parallel Logical. If `TRUE`, parallelizes the grid search
#'   and optimization attempts. See "Parallel Execution" in `@details`.
#' @param n_cores Integer. The number of parallel cores to use (>= 1). If `NULL` (default),
#'   uses `future::availableCores() - 1`.
#'
#' @return A `list` containing the fitted parameters (`beta`), fit statistics
#'   (`negLogLik`, `BIC`, `AICc`, `AIC`), number of parameters (`k`),
#'   number of trials (`N`), and information on fixed parameters.
#'
#' @examples
#' # 1. Load the example data
#' \dontrun{
#' data(SATdata) # Assuming SATdata is an exported dataset
#'
#' # 2. Use a subset of the data for a single fit
#' # (as if it's one participant)
#' single_subject_data <- head(SATdata, 1000)
#' single_subject_data$condition <- factor(single_subject_data$condition)
#'
#' # 3. Define manipulations and fixed parameters
#' manipulations <- list(v ~ condition, a ~ SAT)
#' fixed_params <- list(z = 0.5, svis = 1, d = 0, s = 1)
#'
#' # 4. Fit the model
#' # (Using low attempts/restarts for a quick example)
#' # 'n_cores = 4' will parallelize the n_attempts
#' fit <- fit_rtconf_formula(
#'   single_subject_data,
#'   model = "dynaViTE",
#'   manipulations = manipulations,
#'   fixed = fixed_params,
#'   n_ratings = 6,
#'   sym_thetas = TRUE,
#'   logging = FALSE,
#'   parallel = TRUE,
#'   n_cores = 4,
#'   opts = list(n_attempts = 4, n_restarts = 1)
#' )
#'
#' # 5. View the results
#' print(fit)
#' }
#'
#' @export
#' @md
#' @rdname fit_rtconf_formula
#' @author Sebastian Hellmann, Manuel Methasani
#' @references
#' Hellmann, S., Zehetleitner, M., & Rausch, M. (2023). Simultaneous modeling
#' of choice, confidence, and response time in visual perception.
#' *Psychological Review*. Advance online publication.
#' <https://doi.org/10.1037/rev0000411>
fit_rtconf_formula <- function(
  data,
  model = "dynWEV",
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
  do.call(validate_rtconf_args, args)

  base_context <- args |>
    get_base_context() |>
    setup_main_logging() |>
    setup_parallel()

  if (base_context$parallel) on.exit(plan(sequential), add = TRUE)

  fit_rtconf_formula_dispatcher(base_context)
}

#' @keywords internal
fit_rtconf_formula_dispatcher <- function(base_context) {
  context <- base_context |>
    process_input_data() |>
    get_model_params() |>
    build_model_matrix()

  switch(context$model_type,
    "dynWEV" = fitting_dynwev_formula(context),
    "RM" = stop(sprintf("Model: %s not yet implemented", context$model)),
    stop(sprintf(
      "Model not known. Supported models are: %s",
      paste(DYNWEV_MODELS, collapse = ", ")
    ))
  )
}