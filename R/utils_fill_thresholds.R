#' @keywords internal
fill_thresholds <- function(
  beta,
  used_ratings,
  initial_n_ratings,
  sym_thetas,
  min_conf = -1e+24,
  max_conf = 1e+24
) {
  # store the original parameters that are *not* thresholds
  original_beta <- beta[!grepl("^(d)?theta", names(beta))]

  if (sym_thetas) {
    # case 1: symmetric thresholds (only one set to fill)
    thetas <- fill_single_theta_set(beta, "theta", used_ratings, initial_n_ratings, min_conf, max_conf)
    new_beta <- c(original_beta, thetas)
  } else {
    # case 2: asymmetric thresholds (must fill both 'Upper' and 'Lower')
    thetas_upper <- fill_single_theta_set(beta, "thetaUpper", used_ratings, initial_n_ratings, min_conf, max_conf)
    thetas_lower <- fill_single_theta_set(beta, "thetaLower", used_ratings, initial_n_ratings, min_conf, max_conf)
    new_beta <- c(original_beta, thetas_upper, thetas_lower)
  }

  new_beta
}

#' @keywords internal
#' @noRd
fill_single_theta_set <- function(
  beta,
  theta_prefix,
  used_ratings,
  initial_n_ratings,
  min_conf,
  max_conf
) {
  n_thetas <- initial_n_ratings - 1
  full_thetas <- setNames(rep(NA_real_, n_thetas), paste0(theta_prefix, 1:n_thetas))

  fitted_theta_names <- grep(pattern = paste0("^", theta_prefix, "[0-9]"), names(beta), value = TRUE)

  # get the indices of the thresholds that were actually used
  # (e.g., if ratings 1, 3, 4 were used, this is c(1, 3))
  used_theta_indices <- used_ratings[used_ratings < initial_n_ratings]
  used_theta_names <- paste0(theta_prefix, used_theta_indices)

  # fill in the NAs with the fitted values
  full_thetas[used_theta_names] <- beta[fitted_theta_names]

  min_used <- min(used_ratings)
  if (min_used > 1) {
    # if the subject never used rating 1 (or 1 and 2),
    # fill the lowest thresholds with the minimum value
    low_indices <- seq_len(min_used - 1)
    full_thetas[low_indices] <- min_conf
  }

  max_used <- max(used_ratings)
  if (max_used < initial_n_ratings && max_used <= n_thetas) {
    # if the subject never used the highest rating,
    # fill the highest thresholds with the maximum value
    high_indices <- max_used:n_thetas
    full_thetas[high_indices] <- max_conf
  }

  # fill any remaining internal NAs (e.g., if rating 2 was missed but 1 and 3 were used)
  # by carrying forward the last valid threshold
  for (i in seq_along(full_thetas)) {
    if (i > 1 && is.na(full_thetas[i])) {
      full_thetas[i] <- full_thetas[i - 1]
    }
  }

  full_thetas
}
