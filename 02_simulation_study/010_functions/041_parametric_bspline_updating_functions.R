
# Additional Updating functions for param_bayes


emax_update <- function(y_tilde, h_theta0, emax_sq, emax_mu, sigma_sq){
  emax_sq_new <- ( (emax_sq * t(h_theta0) %*% h_theta0 + sigma_sq) /
                     (sigma_sq * emax_sq))^(-1)
  emax_mu_new <- (emax_sq * t(y_tilde) %*% h_theta0 + emax_mu*sigma_sq) /
    (sigma_sq * emax_sq) * (emax_sq_new)

  emax_new <- emax_mu_new + sqrt(emax_sq_new)*(rnorm(1))

  return(emax_new)
}

#-------------------------------------------------------------------------------

update_ED50_param_bspline <- function(x, ED50_curr, ED50_mu, ED50_v_prop, ED50_v, y_tilde,
                                    h, sigma_sq_tilde){

  n <- length(x)

  lower <- pnorm(0, mean = ED50_curr, sd = sqrt(ED50_v_prop))
  u <- runif(1, lower, 1)
  ED50_new <- qnorm(u, ED50_curr, sd = sqrt(ED50_v_prop))


  h_theta0_ED50 = function(ED50) sigEmax(x, 0, 1, ED50, h)

  y_tilde <- as.vector(y_tilde)

  # log posterior probability of current
  p_curr_log <- dmvnorm(y_tilde, mean = h_theta0_ED50(ED50_curr),
                        sigma = sigma_sq_tilde * diag(n), log = TRUE) +
    # prior:
    dnorm(ED50_curr, mean = ED50_mu, sd = sqrt(ED50_v), log = TRUE)

  # log posterior probability of new
  p_new_log <- dmvnorm(y_tilde, mean = h_theta0_ED50(ED50_new),
                       sigma = sigma_sq_tilde * diag(n), log = TRUE) +
    # prior:
    dnorm(ED50_new, mean = ED50_mu, sd = sqrt(ED50_v), log = TRUE)

  # proposal current given new and new given current
  prop_curr_new_log <- dnorm(ED50_curr, mean = ED50_new, sd = sqrt(ED50_v), log = TRUE)
  prop_new_curr_log <- dnorm(ED50_new, mean = ED50_curr, sd = sqrt(ED50_v), log = TRUE)

  # log Hastings ratio
  HR_log <- p_new_log + prop_curr_new_log - (p_curr_log + prop_new_curr_log)

  # Accept or reject porposal
  if (exp(HR_log)>1 || runif(1) < exp(HR_log)){
    return(ED50_new)
  } else {
    return(ED50_curr)
  }

}

#-------------------------------------------------------------------------------

update_h_param_bspline <- function(x, h_curr, ED50, h_mu, h_v, h_v_prop, y_tilde, sigma_sq_tilde){

  n <- length(x)
  h_new <- rlnorm(1, mean = log(h_curr) - h_v/2, sdlog = sqrt(h_v_prop))

  h_theta0_h = function(h) sigEmax(x, 0, 1, ED50, h)

  y_tilde <- as.vector(y_tilde)

  # log posterior probability of current
  p_curr_log <- dmvnorm(y_tilde, mean = h_theta0_h(h_curr),
                        sigma = sigma_sq_tilde * diag(n), log = TRUE) +
    # prior:
    dlnorm(h_curr, mean = h_mu, sd = sqrt(h_v), log = TRUE)

  # log posterior probability of new
  p_new_log <- dmvnorm(y_tilde, mean = h_theta0_h(h_new),
                       sigma = sigma_sq_tilde * diag(n), log = TRUE) +
    # prior:
    dlnorm(h_new, mean = h_mu, sd = sqrt(h_v), log = TRUE)

  # proposal current given new and new given current
  prop_curr_new_log <- dlnorm(h_curr, mean = log(h_new) - h_v/2, sd = sqrt(h_v), log = TRUE)
  prop_new_curr_log <- dlnorm(h_new, mean = log(h_curr) - h_v/2, sd = sqrt(h_v), log = TRUE)

  # log Hastings ratio
  HR_log <- p_new_log + prop_curr_new_log - (p_curr_log + prop_new_curr_log)

  # Accept or reject proposal
  if (exp(HR_log)>1 || runif(1) < exp(HR_log)){
    return(h_new)
  } else {
    return(h_curr)
  }

}

#-------------------------------------------------------------------------------

update_d_param_bspline <- function(x, d_curr, d_mu, d_v_prop, d_v,
                       y_tilde, sigma_sq_tilde){

  n <- length(x)

  d_new <- rnorm(1, d_curr, sd = sqrt(d_v_prop))

  power_model <- function(x, a, b, d) a + b*x^d

  h_theta0_d = function(d) power_model(x, 0, 1, d)
  h_theta0_d <- Vectorize(h_theta0_d)

  y_tilde <- as.vector(y_tilde)

  # log posterior probability of current
  p_curr_log <- dmvnorm(y_tilde, mean = h_theta0_d(d_curr),
                        sigma = sigma_sq_tilde * diag(n), log = TRUE) +
    # prior:
    dnorm(d_curr, mean = d_mu, sd = sqrt(d_v), log = TRUE)

  # log posterior probability of new
  p_new_log <- dmvnorm(y_tilde, mean = h_theta0_d(d_new),
                       sigma = sigma_sq_tilde * diag(n), log = TRUE) +
    # prior:
    dnorm(d_new, mean = d_mu, sd = sqrt(d_v), log = TRUE)

  # proposal current given new and new given current
  prop_curr_new_log <- dnorm(d_curr, mean = d_new, sd = sqrt(d_v_prop), log = TRUE)
  prop_new_curr_log <- dnorm(d_new, mean = d_curr, sd = sqrt(d_v_prop), log = TRUE)

  # log Hastings ratio
  HR_log <- p_new_log + prop_curr_new_log - (p_curr_log + prop_new_curr_log)

  # Accept or reject porposal
  if (exp(HR_log)>1 || runif(1) < exp(HR_log)){
    return(d_new)
  } else {
    return(d_curr)
  }

}

