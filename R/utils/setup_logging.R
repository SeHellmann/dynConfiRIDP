#' @importFrom logger
#'   layout_glue_generator
#'   log_layout log_appender
#'   appender_file
#'   log_threshold
#'   log_info
#'   log_success
#'   DEBUG
#' @keywords internal
#' @noRd
setup_logging <- function(context) {
  if (context$logging) {
    participant <- 999
    cols <- names(context$data)
    if ("sbj" %in% cols) {
      sbjcol <- "sbj"
    } else if ("participant" %in% cols) {
      sbjcol <- "participant"
    } else {
      sbjcol <- NULL
    }

    if (!is.null(sbjcol)) {
      unique_participants <- unique(context$data[[sbjcol]])
      if (length(unique_participants) == 1) {
        participant <- unique_participants[1]
      }
    }

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
    logdir <- file.path("autosave", paste0("fit", context$model))
    dir.create(logdir, showWarnings = FALSE)

    logfile <- file.path(logdir, paste0("logging_", context$model, ".txt"))

    log_layout <- logger::layout_glue_generator(
      format = paste(
        "{level} [{time}] on process {pid} {fn} for participant",
        participant,
        "and model",
        context$model,
        ": {msg}"
      )
    )

    logger::log_layout(log_layout)
    logger::log_appender(logger::appender_file(file = logfile), index = 2)
    logger::log_threshold(logger::DEBUG, index = 2)
    logger::log_threshold(logger::DEBUG)

    context$logfile <- logfile

    context
  }
}
