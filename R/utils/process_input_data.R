#' @keywords internal
#' @noRd
process_input_data <- function(
  data,
  fixed,
  sym_thetas,
  n_ratings
) {
  cols <- names(data)

  has_rt <- "rt" %in% cols
  has_rating <- "rating" %in% cols

  has_response <- "response" %in% cols
  has_correct <- "correct" %in% cols
  has_stimulus <- "stimulus" %in% cols

  stopifnot(has_rt, "`data` must contain an `rt` column")
  stopifnot(has_rating, "`data` must contain a `rating` column")
  # either 2 of stimulus - response - correct
  # if all 3 => ignore correct
  # stimulus c(-1, 1)
  # response c("lower", "upper")
  # correct c(0, 1)
  # => better validation, clear dependent_vars data.frame spec
  stopifnot(has_response || has_correct, "`data` must contain at least `response` or `correct`")

  rt <- data[["rt"]]
  rating <- data[["rating"]]
  stimulus <- if (has_stimulus) data[["stimulus"]] else NULL

  if (has_response) {
    response <- data[["response"]]
  } else if (has_correct) {
    correct <- data[["correct"]]
    if (has_stimulus) {
      stim_levels <- sort(unique(stimulus))
      if (length(stim_levels) != 2) {
        stop(sprintf(
          "`stimulus` must have exactly 2 unique values, found: %s",
          paste(stim_levels, collapse = ", ")
        ))
      }
      stimulus <- ifelse(stimulus == stim_levels[1], -1, 1)
      response <- ifelse(stimulus * (-1)^correct == 1, -1, 1)
    } else {
      is_z_fixed <- ("z" %in% names(fixed)) && (fixed[["z"]] == 0.5)
      has_equal_bounds <- ("a" %in% names(fixed) && "b" %in% names(fixed) && (fixed[["a"]] == fixed[["b"]]))
      if (!is_z_fixed && sym_thetas && !has_equal_bounds) {
        warning(
          "Only `correct` provided, no `stimulus`: bias cannot be estimated\n",
          "Recommended: fix `z = 0.5`, `sym_thetas = TRUE`, or make sure that `a` and `b` are equal"
        )
      }
      response <- correct
    }
  }

  # pretty sure it should be response == 1 here
  response <- if (all(response %in% c(0, 1))) as.logical(response) else response == 1

  if (!is.integer(rating)) rating <- as.integer(as.factor(rating))

  if (is.null(n_ratings)) n_ratings <- max(rating)

  if (n_ratings < 2) {
    stop(sprintf("There must be at least two unique rating levels\n`n_ratings`: %s", n_ratings))
  }

  if (any(rating == 0) && max(length(unique(rating)), max(rating) + 1 - min(rating)) == n_ratings) {
    rating <- rating + 1
  }

  if (length(unique(rating)) < n_ratings) {
    rating <- as.integer(as.factor(rating))
    n_ratings <- length(unique(rating))
  }

  dependent_vars <- data.frame(rating = rating, response = response, rt = rt)
  if (has_stimulus) dependent_vars$stimulus <- stimulus

  list(dependent_vars = dependent_vars, n_ratings = n_ratings)
}

