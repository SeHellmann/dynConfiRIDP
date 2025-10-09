#' @export
fit_rtconf_formula <- function(
  data,
  model = "dynWEV",
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
  do.call(validate_rtconf_args, args)

  context <- args |>
    fit_context() |>
    setup_logging() |>
    setup_parallel() |>
    get_model_params() |>
    process_input_data() |>
    build_model_matrix()

  if (context$parallel) on.exit(future::plan(sequential), add = TRUE)

  return(switch(context$model_type,
    "dynWEV" = fitting_dynwev_formula(context),
    "RM" = stop(sprintf("Model: %s not yet implemented", context$model)),
    stop(sprintf(
      "Model not known. Supported models are: %s",
      paste(DYNWEV_MODELS, collapse = ", ")
    ))
  ))
}
