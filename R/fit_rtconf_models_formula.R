#' @importFrom parallel
#'   detectCores
#'   makeCluster
#'   stopCluster
#'   clusterExport
#'   clusterApplyLB
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
  parallel_mode = "none",
  n_cores = NULL
) {
  validate_rtconf_models_args(
    models,
    optim_method,
    manipulations,
    parallel_mode,
    n_cores
  )

  ### Determine number of jobs, i.e. model-participant-combinations
  sbjcol <- c("subject", "participant", "sbj")[which(c("subject", "participant", "sbj") %in% names(data))]
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
  n_jobs <- length(models) * length(subjects)

  parallel_subject <- FALSE
  parallel_model <- FALSE
  n_cores_subject <- NULL
  n_cores_model <- NULL
  switch(parallel_mode,
    "none" = {},
    "subject" = {
      parallel_subject <- TRUE
      n_cores_subject <- if (is.null(n_cores)) detectCores() - 1 else n_cores
    },
    "model" = {
      parallel_model <- TRUE
      n_cores_model <- if (is.null(n_cores)) min(detectCores() - 1, n_jobs) else n_cores
    },
    "both" = {
      parallel_subject <- TRUE
      parallel_model <- TRUE

      if (is.null(n_cores)) {
        n_cores_subject <- detectCores() - 1
        n_cores_model <- min(detectCores() - 1, n_jobs)
      } else {
        n_cores_subject <- n_cores[1]
        n_cores_model <- n_cores[2]
      }
    }
  )

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
      logging = logging,
      parallel = parallel_subject,
      n_cores = n_cores_subject
    )
    # not completely sure about the return yet
    res$model <- job$model
    res$sbj <- job$sbj

    res
  }

  if (parallel_model) {
    clmodels <- makeCluster(type = "SOCK", n_cores_model)
    clusterExport(clmodels, c(
      "data",
      "optim_method",
      "fixed",
      "n_ratings",
      "restr_tau",
      "sym_thetas",
      "precision",
      "opts",
      "grid_search",
      "logging",
      "parallel_subject",
      "n_cores_subject",
      "call_fitfct"
    ), envir = environment())
    on.exit(try(stopCluster(clmodels), silent = TRUE))
    res <- clusterApplyLB(clmodels, jobs_list, fun = call_fitfct)
    stopCluster(clmodels)
  } else {
    res <- lapply(jobs_list, call_fitfct)
  }

  res
}
