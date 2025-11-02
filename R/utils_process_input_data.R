#' @keywords internal
process_input_data <- function(context) {
  cols <- names(context$data)

  rt <- context$data[["rt"]]
  rating <- context$data[["rating"]]

  response <- context$data[["response"]]
  stimulus <- context$data[["stimulus"]]
  correct <- context$data[["correct"]]

  validate_input_data(
    cols,
    rt,
    rating,
    response,
    stimulus,
    correct
  )

  # parse stimulus into -1 / 1 coding
  if (!is.null(stimulus)) {
    stimulus_levels <- sort(unique(stimulus))
    stimulus <- ifelse(stimulus == stimulus_levels[1], -1, 1)
  }

  # parse response into 0 / 1 coding
  # if response isn't provided, infer it from correct and stimulus
  if (!is.null(response)) {
    response_levels <- sort(unique(response))
    if (!all(response_levels == stimulus_levels)) {
      stop(sprintf(
        "`response` and `stimulus` must have the same unique values\nresponse: %s\nstimulus: %s",
        paste(response_levels, collapse = ", "), paste(stimulus_levels, collapse = ", ")
      ))
    }
    response <- ifelse(response == response_levels[1], 0, 1)
  } else if (!is.null(correct)) {
    if (!is.null(stimulus)) {
      # infer response from stimulus and correct
      # response 0 = stimulus -1
      # response 1 = stimulus 1
      response <- ifelse((stimulus * (-1)^correct == 1) + 1, 0, 1)
    } else {
      # if only correct is given, assume 0 = response 0, 1 = response 1
      response <- ifelse(correct == 0, 0, 1)

      # warn user if bias parameters (z, thetas) might not be identifiable
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

  # ensure ratings are 1-indexed
  if (any(rating == 0)) rating <- rating + 1

  all_ratings <- min(rating):max(rating)
  total_ratings <- length(all_ratings)
  used_n_ratings <- length(unique(rating))

  if (is.null(context$n_ratings)) context$n_ratings <- total_ratings

  # handle cases where a subject did not use all possible ratings
  # we fit only the *used* ratings and fill the gaps later
  if (used_n_ratings < context$n_ratings) {
    context$used_ratings <- sort(unique(rating))
    context$initial_n_ratings <- context$n_ratings
    # re-index ratings to be contiguous (e.g., 1, 3, 4 -> 1, 2, 3)
    rating <- as.integer(as.factor(rating))
    context$n_ratings <- used_n_ratings
  }

  context$dependent_vars <- cbind(rt = rt, rating = rating, response = response)
  if (!is.null(stimulus)) context$dependent_vars <- cbind(context$dependent_vars, stimulus = stimulus)

  # maxt0 is the minimum possible non-decision time (min RT)
  context$maxt0 <- min(rt)

  context
}
