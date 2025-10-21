#' @keywords internal
setup_main_logging <- function(context) {
  if (!context$logging) return(context)

  base_dir <- file.path(tempdir(), "rtconf_fits")
  identifier <- "main"
  log_tag <- "[MAIN]:"

  apply_log_settings(context, base_dir, identifier, log_tag)
}

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

  base_dir <- file.path(tempdir(), "rtconf_fits", paste0("fit_", context$model))
  if (participant_id != 999)
    base_dir <- file.path(base_dir, sprintf("subj_%s", participant_id))
  identifier <- paste(context$model, participant_id, sep = "_")
  log_tag <- paste0("[SUBJ: ", participant_id, "]:")

  apply_log_settings(context, base_dir, identifier, log_tag)
}

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

  log_layout(layout_glue_generator(
    format = paste(log_parts[log_parts != ""], collapse = " ")
  ))
  log_appender(appender_file(log_file), index = 2)
  log_threshold(DEBUG, index = 2)

  context$data_file_prefix <- data_file_prefix
  log_info(sprintf("Logging initialized, logs will be saved in: %s", base_dir))

  context
}

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