#' @keywords internal
setup_logging <- function(context) {
  if (!context$logging) return(context)

  participant_id <- extract_participant_id(context$data)
  model_dir <- setup_log_dirs(context$model)
  paths <- get_log_paths(model_dir, context$model, participant_id)

  logger::log_layout(setup_log_layout(context$model, participant_id))
  logger::log_appender(logger::appender_file(paths$log_file), index = 2)
  logger::log_threshold(logger::DEBUG, index = 2)

  context$log_file <- paths$log_file
  context$data_file <- paths$data_file

  logger::log_info(sprintf("Logging initialized: %s", paths$log_file))
  logger::log_info(sprintf("Data file location: %s", paths$data_file))

  context
}

#' @keywords internal
extract_participant_id <- function(data) {
  participant <- 999
  potential_sbj_cols <- c("sbj", "participant", "subject")
  sbjcol <- potential_sbj_cols[potential_sbj_cols %in% names(data)][1]

  if (!is.null(sbjcol)) {
    unique_participants <- unique(data[[sbjcol]])
    if (length(unique_participants) == 1) {
      participant <- unique_participants[1]
    }
  }

  participant
}

#' @keywords internal
setup_log_dirs <- function(model_name) {
  base_dir <- file.path(tempdir(), "rtconf_fits")
  model_dir <- file.path(base_dir, paste0("fit_", model_name))
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  model_dir
}

#' @keywords internal
get_log_paths <- function(model_dir, model_name, participant_id) {
  unique_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  list(
    log_file = file.path(model_dir, sprintf("log_%s_%d_%s.txt", model_name, participant_id, unique_id)),
    data_file = file.path(model_dir, sprintf("data_%s_%d_%s.RData", model_name, participant_id, unique_id))
  )
}

#' @keywords internal
setup_log_layout <- function(model_name, participant_id) {
  logger::layout_glue_generator(
    format = paste(
      "[{format(time, \"%H:%M:%S\")}]",
      "[PID: {pid}]",
      "[{level}]",
      "[{fn}]",
      paste0("[MODEL: ", model_name, "]"),
      paste0("[SUBJ: ", participant_id, "]:"),
      "{msg}"
    )
  )
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
      context$data_file
    )
  )
}
