#' dynConfiRIDP: Dynamic Models for Confidence and Response Time Distributions
#'
#' Provides response time and confidence distributions (density/PDF) for
#' the following models: dynaViTE, dynWEV, 2DSD, IRM and PCRM.
#' The package also provides a formula-based interface for fitting
#' these models to data.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib dynConfiRIDP, .registration = TRUE
#' @import nloptr
#' @import RcppArmadillo
#' @importFrom  Rcpp
#'   sourceCpp
#' @importFrom assertthat
#'   assert_that
#'   is.string
#' @importFrom logger
#'   appender_file
#'   DEBUG
#'   layout_glue_generator
#'   log_appender
#'   log_layout
#'   log_info
#'   log_success
#'   log_warn
#'   log_error
#'   log_threshold
#' @importFrom future
#'   plan
#'   sequential
#'   multisession
#'   multicore
#'   availableCores
#'   supportsMulticore
#'   nbrOfWorkers
#' @importFrom future.apply
#'   future_lapply
#' @importFrom stats
#'   model.matrix
#'   rnorm
#'   setNames
#' @importFrom utils
#'   capture.output
#'   str
## usethis namespace: end
NULL