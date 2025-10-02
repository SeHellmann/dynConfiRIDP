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

  # parse stimulus
  if (!is.null(stimulus)) {
    stimulus_levels <- sort(unique(stimulus))
    stimulus <- ifelse(stimulus == stimulus_levels[1], -1, 1)
  }

  # infer response if needed and parse
  # if no response then both stimulus and correct or only correct (with warning)
  if (!is.null(response)) {
    response_levels <- sort(unique(response))
    response <- ifelse(response == response_levels[1], -1, 1)
  } else if (!is.null(correct)) {
    if (!is.null(stimulus)) {
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

  context$dependent_vars <- cbind(rt = rt, rating = rating, response = response)
  if (!is.null(stimulus)) context$dependent_vars <- cbind(context$dependent_vars, stimulus = stimulus)

  context$maxt0 <- min(rt)

  # thetas_parnames
  base_names <- if (context$sym_thetas) "theta" else c("thetaLower", "thetaUpper")
  if (context$n_ratings <= 2) {
    thetas_parnames <- base_names
  } else {
    first_thetas <- paste0(base_names, "1")
    dtheta_indices <- 2:(context$n_ratings - 1)
    dthetas <- if (context$sym_thetas) {
      paste0("dtheta", dtheta_indices)
    } else {
      paste0(
        "dtheta",
        rep(c("Lower", "Upper"), times = context$n_ratings - 2),
        rep(dtheta_indices, each = 2)
      )
    }
    thetas_parnames <- c(first_thetas, dthetas)
  }
  context$thetas_parnames <- thetas_parnames

  context
}
