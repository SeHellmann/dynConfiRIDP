#' @keywords internal
fill_thresholds <- function(
  beta,
  used_ratings,
  initial_n_ratings,
  sym_thetas,
  min_conf = -1e+24,
  max_conf = 1e+24
) {
  original_beta <- beta[!grepl("^(d)?theta", names(beta))]

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
  full_thetas <- setNames(rep(NA_real_, n_thetas), paste0(theta_prefix, 1:n_thetas))

  fitted_theta_names <- grep(pattern = paste0("^", theta_prefix, "[0-9]"), names(beta), value = TRUE)
  used_theta_indices <- used_ratings[used_ratings < initial_n_ratings]
  used_theta_names <- paste0(theta_prefix, used_theta_indices)
  full_thetas[used_theta_names] <- beta[fitted_theta_names]

  min_used <- min(used_ratings)
  if (min_used > 1) {
    low_indices <- seq_len(min_used - 1)
    full_thetas[low_indices] <- min_conf
  }

  max_used <- max(used_ratings)
  if (max_used < initial_n_ratings && max_used <= n_thetas) {
    high_indices <- max_used:n_thetas
    full_thetas[high_indices] <- max_conf
  }

  for (i in seq_along(full_thetas)) {
    if (i > 1 && is.na(full_thetas[i])) {
      full_thetas[i] <- full_thetas[i - 1]
    }
  }

  full_thetas
}
