# sets up the primary logger for the main R process
#' @keywords internal
setup_main_logging <- function(context) {
  if (!context$logging) return(context)

  base_dir <- here::here(".logs", "rtconf_fits")
  identifier <- "main"
  log_tag <- "[MAIN]:"

  apply_log_settings(context, base_dir, identifier, log_tag)
}

# sets up a subject-specific logger, creating a dedicated
# log file for a single participant-model fit
#' @keywords internal
setup_subject_logging <- function(context) {
  if (!context$logging) return(context)

  participant_id <- {
    cols <- c("sbj", "participant", "subject")
    sbjcol <- cols[cols %in% names(context$data)][1]

    if (is.na(sbjcol)) {
      999
    } else {
      unique_participants <- unique(context$data[[sbjcol]])
      if (length(unique_participants) == 1) unique_participants[1] else 999
    }
  }

  base_dir <- here::here(".logs", "rtconf_fits", paste0("fit_", context$model))
  if (participant_id != 999) {
    base_dir <- file.path(base_dir, sprintf("subj_%s", participant_id))
  }

  identifier <- paste(context$model, participant_id, sep = "_")
  log_tag <- paste0("[SUBJ: ", participant_id, "]:")

  apply_log_settings(context, base_dir, identifier, log_tag)
}

# used to apply logging settings (layout, file) to a
# parallel worker (e.g., in a future_lapply call)
# this ensures that parallel processes also write to the log file
#' @keywords internal
setup_worker_logging <- function(log_config) {
  if (is.null(log_config)) return(invisible())

  log_layout(log_config$layout, index = 2)
  log_appender(appender_file(log_config$log_file), index = 2)
  log_threshold(log_config$threshold, index = 2)
}

# helper function that creates the log directory and defines
# the log message format (layout)
#' @keywords internal
apply_log_settings <- function(context, base_dir, identifier, log_tag) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(base_dir, sprintf("log_%s_%s.txt", identifier, timestamp))
  data_file_prefix <- file.path(base_dir, sprintf("data_%s_%s.RData", identifier, timestamp))

  model <- if (!is.null(context$model)) {
    context$model
  } else if (!is.null(context$models) && length(context$models) == 1) {
    context$models[1]
  } else {
    NULL
  }

  model_tag <- if (!is.null(model)) {
    paste0("[MODEL: ", model, "]")
  } else {
    ""
  }

  log_parts <- c(
    "[{format(time, \"%H:%M:%S\")}]",
    "[PID: {pid}]",
    "[{level}]",
    "[{fn}]",
    model_tag,
    log_tag,
    "{msg}"
  )
  layout <- layout_glue_generator(
    format = paste(log_parts[log_parts != ""], collapse = " ")
  )

  log_layout(layout, index = 1)
  log_layout(layout, index = 2)
  log_appender(appender_file(log_file), index = 2)
  log_threshold(DEBUG, index = 2)

  # store the log config in the context so parallel workers can access it
  context$log_config <- list(
    layout = layout,
    log_file = log_file,
    threshold = DEBUG
  )
  context$data_file_prefix <- data_file_prefix

  log_info(sprintf("Logging initialized, logs will be saved in: %s", base_dir))

  context
}

# saves intermediate or final optimization states (.RData files)
# to the .logs/ directory for debugging
#' @keywords internal
save_optimization_state <- function(context, state_data, type = c("grid_search", "optimization", "final")) {
  if (!context$logging) return(invisible())

  type <- match.arg(type)
  state_name <- switch(type,
    "grid_search" = "grid_search_state",
    "optimization" = "optim_state",
    "final" = "final_state"
  )
  save(
    state_data,
    file = sub(
      "\\.RData$",
      paste0("_", state_name, ".RData"),
      context$data_file_prefix
    )
  )
}