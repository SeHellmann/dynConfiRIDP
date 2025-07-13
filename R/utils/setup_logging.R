#' @keywords internal
#' @noRd
setup_logging <- function(logging, participant, model) {
  if (!requireNamespace("logger", quietly = TRUE)) {
    warning(
      "Package 'logger' is not installed but needed for logging.\n",
      "Process continues without logging.\n",
      "Interrupt and install 'logger' if logging is needed.",
      immediate. = TRUE
    )
    return(invisible(NULL))
  }

  dir.create("autosave", showWarnings = FALSE)
  logdir <- file.path("autosave", paste0("fit", model))
  dir.create(logdir, showWarnings = FALSE)

  logfile <- file.path(logdir, paste0("logging_", model, ".txt"))

  log_layout <- logger::layout_glue_generator(
    format = paste(
      "{level} [{time}] on process {pid} {fn} for participant",
      participant,
      "and model",
      model,
      ": {msg}"
    )
  )

  logger::log_layout(log_layout)
  logger::log_appender(logger::appender_file(file = logfile), index = 2)
  logger::log_threshold(logger::DEBUG, index = 2)
  logger::log_threshold(logger::DEBUG)

  logfile
}
