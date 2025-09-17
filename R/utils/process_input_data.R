#' @keywords internal
#' @noRd
process_input_data <- function(context) {
  cols <- names(context$data)

  has_rt <- "rt" %in% cols
  has_rating <- "rating" %in% cols

  has_response <- "response" %in% cols
  has_stimulus <- "stimulus" %in% cols
  has_correct <- "correct" %in% cols

  stopifnot(has_rt, "`data` must contain a `rt` column")
  stopifnot(has_rating, "`data` must contain a `rating` column")
  # if no response then both stimulus and correct or only correct (with warning)
  stopifnot(has_response || has_correct, "`data` must contain at least `response` or `correct`" )

  rt <- context$data[["rt"]]
  rating <- context$data[["rating"]]

  response <- if (has_response) context$data[["response"]] else NULL
  stimulus <- if (has_stimulus) context$data[["stimulus"]] else NULL
  correct <- if (has_correct) context$data[["correct"]] else NULL

  validate_input_data(
    rt,
    rating,
    response,
    stimulus,
    correct
  )

  # parse stimulus
  if (has_stimulus) {
    stimulus_levels <- sort(unique(stimulus))
    stimulus <- ifelse(stimulus == stimulus_levels[1], -1, 1)
  }

  # infer response if needed and parse
  if (has_response) {
    response_levels <- sort(unique(response))
    response <- ifelse(response == response_levels[1], -1, 1)
  } else if (has_correct) {
    if (has_stimulus) {
      response <- ifelse(stimulus * (-1)^correct == 1, -1, 1)
    } else {
      # get response from correct
      response <- ifelse(correct == 0, -1, 1)

      is_z_fixed <- ("z" %in% names(context$fixed)) && (context$fixed[["z"]] == 0.5)
      has_equal_bounds <- (
        "a" %in% names(context$fixed) &&
        "b" %in% names(context$fixed) &&
        (context$fixed[["a"]] == context$fixed[["b"]])
      )
      if (!is_z_fixed && context$sym_thetas && !has_equal_bounds) {
        warning(
          "Only `correct` provided, no `stimulus`: bias cannot be estimated\n",
          "Recommended: fix `z = 0.5`, `sym_thetas = TRUE`, or make sure that `a` and `b` are equal"
        )
      }
    }
  }

  if (any(rating == 0)) rating <- rating + 1

  all_ratings <- min(rating):max(rating)
  total_ratings <- length(all_ratings)
  used_n_ratings <- length(unique(rating))

  if (is.null(context$n_ratings)) context$n_ratings <- total_ratings

  if (used_n_ratings < context$n_ratings) {
    context$used_ratings <- sort(unique(rating))
    context$initial_n_ratings <- context$n_ratings
    rating <- as.integer(as.factor(rating))
    context$n_ratings <- used_n_ratings
  }

  context$dependent_vars <- data.frame(rt = rt, rating = rating, response = response)
  if (has_stimulus) context$dependent_vars$stimulus <- stimulus

  context$maxt0 <- min(rt)

  context
}

