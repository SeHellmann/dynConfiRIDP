#' @importFrom logger
#'   layout_glue_generator
#'   log_layout log_appender
#'   appender_file
#'   log_threshold
#'   log_info
#'   log_success
#'   DEBUG
#' @keywords internal
setup_logging <- function(context) {
  if (context$logging) {
    participant <- 999
    potential_sbj_cols <- c("sbj", "participation", "subject")
    sbjcol <- potential_sbj_cols[match(TRUE, potential_sbj_cols %in% names(context$data))]
    if (is.na(sbjcol)) sbjcol <- NULL

    if (!is.null(sbjcol)) {
      unique_participants <- unique(context$data[[sbjcol]])
      if (length(unique_participants) == 1) {
        participant <- unique_participants[1]
      }
    }

    logdir <- file.path(tempdir(), paste0("rtconf_fits"))
    dir.create(logdir, showWarnings = FALSE)

    unique_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
    logfile <- file.path(logdir, paste0("log_", context$model, "_", unique_id, ".txt"))

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
    logger::log_info(sprintf("Logging to file: %s", logfile))

    context$logfile <- logfile
  }

  context
}
