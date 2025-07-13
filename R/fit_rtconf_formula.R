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
  validate_rtconf_args(
    data,
    model,
    optim_method,
    fixed,
    manipulations,
    n_ratings,
    restr_tau,
    sym_thetas,
    precision,
    opts,
    grid_search,
    logging,
    parallel,
    n_cores
  )

  #### setup logging ####
  if (logging) {
    participant <- 999
    cols <- names(data)
    if ("sbj" %in% cols) {
      sbjcol <- "sbj"
    } else if ("participant" %in% cols) {
      sbjcol <- "participant"
    } else {
      sbjcol <- NULL
    }

    if (!is.null(sbjcol)) {
      unique_participants <- unique(data[[sbjcol]])
      if (length(unique_participants) == 1) {
        participant <- unique_participants[1]
      }
    }
    logfile <- setup_logging(logging, model, participant)
  }

  #### get all params for the given model ####
  model_params <- get_model_params(model, fixed, manipulations)
  # maybe also unnecessary do destructure the list like this
  model_type <- model_params$model_type
  parnames <- model_params$parnames
  fixed <- model_params$fixed
  fixed_parnames <- model_params$fixed_parnames
  manipulations <- model_params$manipulations
  manipulated_parnames <- model_params$manipulated_parnames
  const_parnames <- model_params$const_parnames

  #### fill in default optimizer options if missing ####
  missing_opts <- !(names(DEFAULT_OPTS) %in% names(opts))
  opts <- c(opts, DEFAULT_OPTS[missing_opts])

  #### process data input ####
  processed_input_data <- process_input_data(data, fixed, sym_thetas, n_ratings)
  dependent_vars <- processed_input_data$dependent_vars
  n_ratings <- processed_input_data$n_ratings

  #### build model_matrix ####
  model_matrix_result <- build_model_matrix(data, manipulated_parnames, manipulations)
  model_matrix <- model_matrix_result$model_matrix
  fit_params_cols <- model_matrix_result$fit_params_cols
  fit_beta_parnames <- model_matrix_result$fit_beta_parnames
  fit_parnames <- model_matrix_result$fit_parnames

  return(switch(model_type,
    "dynWEV" = fitting_dynwev_formula(),
    "RM" = fitting_rm_formula(),
    stop(sprintf(
      "Model not known. Model must be one of: %s",
      paste(MODELS, collapse = ", ")
    ))
  ))
}
