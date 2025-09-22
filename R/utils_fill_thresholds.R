#' @keywords internal
fill_thresholds <- function(
  beta,
  used_ratings,
  initial_n_ratings,
  sym_thetas,
  min_conf = -1e+24,
  max_conf = 1e+24
) {
  original_beta <- beta[!grepl("^theta", names(beta))]

  if (sym_thetas) {
    thetas <- fill_single_theta_set(beta, "theta", used_ratings, initial_n_ratings, min_conf, max_conf)
    new_beta <- c(original_beta, thetas)
  } else {
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
  full_thetas <- rep(NA, n_thetas)
  names(full_thetas) <- paste0(theta_prefix, 1:n_thetas)

  fitted_theta_names <- grep(pattern = paste0("^", theta_prefix, "[0-9]"), names(beta), value = TRUE)
  used_theta_names <- paste0(theta_prefix, used_ratings[used_ratings < initial_n_ratings])
  full_thetas[used_theta_names] <- beta[fitted_theta_names]

  if (min(used_ratings) > 1) {
    unused_low_names <- paste0(theta_prefix, 1:(min(used_ratings) - 1))
    fill_value <- min(min_conf, min(full_thetas, na.rm = TRUE))
    full_thetas[unused_low_names] <- fill_value
  }
  if (max(used_ratings) < initial_n_ratings) {
    unused_high_names <- paste0(theta_prefix, max(used_ratings):n_thetas)
    fill_value <- max(max_conf, max(full_thetas, na.rm = TRUE))
    full_thetas[unused_high_names] <- fill_value
  }

  for (i in seq_along(full_thetas)) {
    if (i > 1 && is.na(full_thetas[i])) {
      full_thetas[i] <- full_thetas[i - 1]
    }
  }

  full_thetas
}
