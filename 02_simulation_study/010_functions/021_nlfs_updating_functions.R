
# beta update ------------------------------------------------------------------

beta_update <- function(y, Sigma_prec, tau_sq, sigma_sq, alpha, X_mat){

  # Posterior Covariance matrix without scaling with 1 / sigma_sq
  Sigma_new <-tryCatch({chol2inv(chol(crossprod(X_mat) + 1/tau_sq *Sigma_prec))},
                       error = function(cond){
                         message("beta update error in : chol(crossprod(X_mat) + 1/tau_sq *Sigma_prec)")
                         message("Add small constant to catch it")

                         return(chol2inv(chol(crossprod(X_mat) + 1/tau_sq *Sigma_prec + 1e-7*diag(nrow(Sigma_prec)))))

                         })#instead of solve try chol2inv(chol())

  # Complete posterior mu
  mu_new <- Sigma_new %*% t(X_mat) %*% (y - alpha)

  # Cholesky decomposition for more efficient drawing of fully cond. post.
  Sigma_new_Ut <- chol(Sigma_new)

  # Now factor sigma (not sigma_sq because Cov^(1/2)) for scaling back:
  beta <- mu_new + sqrt(sigma_sq)*t(Sigma_new_Ut) %*% rnorm(nrow(Sigma_new))

  return(beta)

}

# interc update ----------------------------------------------------------------

interc_update <- function(y, X_mat, n, alpha_mu, alpha_v, sigma_sq, beta){

  alpha_v_new <- (sigma_sq * alpha_v) / (n * alpha_v + sigma_sq)
  alpha_mu_new <- alpha_v_new * ((n * mean(y) -  sum(X_mat %*% beta) ) / sigma_sq +
                                   alpha_mu / alpha_v)
  alpha <- alpha_mu_new + sqrt(alpha_v_new) * rnorm(1)

  return(alpha)

}

# sigma_sq update --------------------------------------------------------------

sigma_sq_update <- function(sigma_sq_a_new, sigma_sq_b, n, k, tau_sq,
         beta, Sigma_prec, alpha, X_mat, y){

  RSS <- sum( (y - (alpha + X_mat %*% beta) )^2 )
  sigma_sq_b_new <- 0.5 * (RSS + tau_sq^(-1) * t(beta) %*% Sigma_prec %*% beta) +
    sigma_sq_b

  sigma_sq <- invgamma::rinvgamma(1, shape = sigma_sq_a_new, rate = sigma_sq_b_new)

  return(sigma_sq)
}

# tau_sq update for shrinakge == "own_slice" ---------------------------------------


tau_sq_update <- function(beta, old_tau_sq, omega_a, omega_b, k, Sigma_prec,
                          sigma_sq, lower_thres = sqrt(0.001), upper_thres = sqrt(10)){
  old_tau <- sqrt(old_tau_sq)

  c1 <- (-k/2 + omega_b - 0.5) # removed: k-d0
  c2 <- (- omega_a - omega_b)
  c3 <- -1/(2*sigma_sq)* t(beta) %*% Sigma_prec %*% beta
  log_post_tau <- function(tau){
    c1 *log(tau^2) +
      c2 *log(1+tau^2) +
      c3 * 1/(tau^2)
  }
  log_post_tau <- Vectorize(log_post_tau)

  # Find f(curr_tau_sq) to draw from a vertical line
  v <- log_post_tau(old_tau)
  # To be safe: back on regular, non-log level to draw the unif for the slice
  v_exp <- exp(v)
  z_exp = runif(1, 0, v_exp)
  # Now back on log-scale as it's log-posterior
  z = log(z_exp)

  # Search for slice endings, i.e. where log_post_tau has function value z

  levelled_g <- function(tau) log_post_tau(tau) - z
  # find maximum of levelled_g
  #g_amax <- optimize(levelled_g, interval = c(lower_thres, upper_thres), maximum = TRUE)$maximum
  upper <- tryCatch(
    {uniroot(levelled_g, interval = c(old_tau,1e6))$root},
    error = function(cond){
      return(upper_thres)
    })

  lower <-  tryCatch(
    {uniroot(levelled_g, interval = c(1e-10, old_tau))$root},
    error = function(cond){
      return(lower_thres)
    })

  if(lower < lower_thres) lower <- lower_thres
  if(upper > upper_thres) upper <- upper_thres

  #abline(v = c(lower, upper), col = "red")

  if(lower < upper){
    new_tau_sq <- (runif(1, lower, upper))^2
  } else {
    new_tau_sq <- old_tau_sq
  }

  return(new_tau_sq)
}


# tau_sq_update for shrinkage == "half_cauchy" ---------------------------------

tau_sq_update_hc <- function(beta, Sigma_prec, xi, sigma_sq){
  a_new <- (length(beta) + 1) / 2
  b_new <- 1 / (2*sigma_sq) * t(beta) %*% Sigma_prec %*% beta + 1/xi

  res <- rinvgamma(1, shape = a_new, rate = b_new)
  return(res)
}


# Current ED50 update ----------------------------------------------------------

# beta = all_beta[, i]
# tau_sq = all_tau_sq[i]
# sigma_sq = all_sigma_sq[i]
# ED50_curr = all_theta[1, i-1]
# h = all_theta[2, i-1]
# x_grid = d$x
# X = B
# XPX_curr = BIPB
# ED50_mu = theta3_mu
# ED50_v = theta3_sq


ED50_update <- function(y, tau_sq, sigma_sq, ED50_curr,
			nlf_deriv, h, d_par, x_grid, X, Sigma_curr, ED50_mu, ED50_v, ED50_v_prop){

  n <- length(x_grid)
  # Step 1: Proposal distribution
  # Propose ED50_new given ED_50 (truncated positive normal)
  lower <- pnorm(0, mean = ED50_curr, sd = sqrt(ED50_v_prop))

  # repeat this until it works
  try_inv <- 1
  while(try_inv){
    u <- runif(1, lower, 1)
    ED50_new <- qnorm(u, ED50_curr, sd = sqrt(ED50_v_prop))

    # Calculate Hastings ratio
    # Projection matrices and precision matrices
    # P_new <- get_P(x_grid = x_grid, ED50 = ED50_new, h = h)

    F_mat <- nlf_deriv(x_grid, na.omit(c(ED50_new, h, d_par)))
    # Use no more:
    # P <- F_mat %*% solve(t(F_mat) %*% F_mat) %*% t(F_mat)
    # Use qr decomposition instead
    L <- t(suppressWarnings(chol(crossprod(F_mat), pivot = TRUE)))
    d0 <- attr(L, "rank")
    Qt <- forwardsolve(L, t(F_mat[, attr(L, "pivot")]), d0)
    P_new <- crossprod(Qt)
    XPX_new <- t(X) %*%(diag(n) - P_new) %*% X
    XPX_new_inv <- tryCatch({chol2inv(chol(XPX_new))},
                            error = function(cond){
                              return(NA)
                            })
    if(is.matrix(XPX_new_inv)){
      try_inv <- 0
    } else {
      message("ED50 update: chol2inv(chol(XPX_new)) threw error, try again")
      message(paste0("ED50_new was ", ED50_new))
    }
  }

  Sigma_new <- sigma_sq * tau_sq * X %*% XPX_new_inv %*% t(X) + sigma_sq * diag(n)

  # For current: Already calculated


  # log posterior probability of current
  p_curr_log <- dmvnorm(y, mean = rep(0, n), sigma = Sigma_curr, log = TRUE) +
    # prior:
    dnorm(ED50_curr, mean = ED50_mu, sd = sqrt(ED50_v), log = TRUE)

  # log posterior probability of new
  p_new_log <- dmvnorm(y, mean = rep(0, n), sigma = Sigma_new, log = TRUE) +
    # prior:
    dnorm(ED50_new, mean = ED50_mu, sd = sqrt(ED50_v), log = TRUE)

  # proposal current given new and new given current (symmetric, could be ignored)
  prop_curr_new_log <- dnorm(ED50_curr, mean = ED50_new, sd = sqrt(ED50_v_prop), log = TRUE)
  prop_new_curr_log <- dnorm(ED50_new, mean = ED50_curr, sd = sqrt(ED50_v_prop), log = TRUE)

  # log Hastings ratio
  HR_log <- p_new_log + prop_curr_new_log - (p_curr_log + prop_new_curr_log)

  # Accept or reject proposal
  if (exp(HR_log)>1 || runif(1) < exp(HR_log)){
    return(ED50_new)
  } else {
    return(ED50_curr)
  }
}


# Current h update -------------------------------------------------------------

h_update <- function(y, tau_sq, sigma_sq, ED50, d_par, h_curr, x_grid, X, Sigma_curr,
			nlf_deriv,
                     h_mu, h_v, h_v_prop){

  n <- length(x_grid)

  # Heuristic for proposal distribution: (heu_fac* h_curr) as sd.
  try_inv <- 1
  while(try_inv){
    h_new <- rlnorm(1, mean = log(h_curr) - h_v_prop/2, sdlog = sqrt(h_v_prop))


    # Calculate Hastings ratio
    # Projection matrices and precision matrices
    F_mat <- nlf_deriv(x_grid, na.omit(c(ED50, h_new, d_par)))
    # Use no more:
    # P <- F_mat %*% solve(t(F_mat) %*% F_mat) %*% t(F_mat)
    # Use qr decomposition instead
    L <- t(suppressWarnings(chol(crossprod(F_mat), pivot = TRUE)))
    d0 <- attr(L, "rank")
    Qt <- forwardsolve(L, t(F_mat[, attr(L, "pivot")]), d0)
    P_new <- crossprod(Qt)
    XPX_new <- t(X) %*%(diag(n) - P_new) %*% X
    XPX_new_inv <- tryCatch({chol2inv(chol(XPX_new))},
                            error = function(cond){
                              return(NA)
                            })
    if(is.matrix(XPX_new_inv)){
      try_inv <- 0
    } else {
      message("h update: chol2inv(chol(XPX_new)) threw error, try again")
      message(paste0("h_new was ", h_new))
    }
  }
  Sigma_new <- sigma_sq * tau_sq * X %*% XPX_new_inv %*% t(X) + sigma_sq * diag(n)

  # Also for the current: Already calculated

  # log posterior probability of current
  p_curr_log <- dmvnorm(y, mean = rep(0, n), sigma = Sigma_curr, log = TRUE) +
    # prior:
    dlnorm(h_curr, mean = h_mu, sd = sqrt(h_v), log = TRUE)

  # log posterior probability of new
  p_new_log <- dmvnorm(y, mean = rep(0, n), sigma = Sigma_new, log = TRUE) +
    # prior:
    dlnorm(h_new, mean = h_mu, sd = sqrt(h_v), log = TRUE)

  # proposal current given new and new given current
  prop_curr_new_log <- dlnorm(h_curr, mean = log(h_new) - h_v/2, sd = sqrt(h_v_prop), log = TRUE)
  prop_new_curr_log <- dlnorm(h_new, mean = log(h_curr) - h_v/2, sd = sqrt(h_v_prop), log = TRUE)

  # log Hastings ratio
  HR_log <- p_new_log + prop_curr_new_log - (p_curr_log + prop_new_curr_log)

  # Accept or reject porposal
  if (exp(HR_log)>1 || runif(1) < exp(HR_log)){
    return(h_new)
  } else {
    return(h_curr)
  }
}


# d_update in power model ------------------------------------------------------

# f(x)=a + b*x^d
# d is only non-linear parameter in power-model
# normal prior is assumed for d

d_update <- function(y, tau_sq, sigma_sq, d_curr, ED50, h, nlf_deriv,
                     x_grid, X, Sigma_curr, d_mu, d_v, d_v_prop){

  n <- length(x_grid)
  # Step 1: Proposal distribution
  # Propose ED50_new given ED_50 (truncated positive normal)

  # repeat this until it works
  try_inv <- 1
  while(try_inv){
    d_new <- rnorm(1, d_curr, sd = sqrt(d_v_prop))

    # Calculate Hastings ratio
    # Projection matrices and precision matrices
    # P_new <- get_P(x_grid = x_grid, ED50 = ED50_new, h = h)

    F_mat <- nlf_deriv(x_grid, na.omit(c(ED50, h, d_new)))
    # Use no more:
    # P <- F_mat %*% solve(t(F_mat) %*% F_mat) %*% t(F_mat)
    # Use qr decomposition instead
    L <- t(suppressWarnings(chol(crossprod(F_mat), pivot = TRUE)))
    d0 <- attr(L, "rank")
    Qt <- forwardsolve(L, t(F_mat[, attr(L, "pivot")]), d0)
    P_new <- crossprod(Qt)
    XPX_new <- t(X) %*%(diag(n) - P_new) %*% X
    XPX_new_inv <- tryCatch({chol2inv(chol(XPX_new))},
                            error = function(cond){
                              return(NA)
                            })
    if(is.matrix(XPX_new_inv)){
      try_inv <- 0
    } else {
      message("d update: chol2inv(chol(XPX_new)) threw error, try again")
      message(paste0("d_new was ", d_new))
    }
  }

  Sigma_new <- sigma_sq * tau_sq * X %*% XPX_new_inv %*% t(X) + sigma_sq * diag(n)

  # For current: Already calculated


  # log posterior probability of current
  p_curr_log <- dmvnorm(y, mean = rep(0, n), sigma = Sigma_curr, log = TRUE) +
    # prior:
    dnorm(d_curr, mean = d_mu, sd = sqrt(d_v), log = TRUE)

  # log posterior probability of new
  p_new_log <- dmvnorm(y, mean = rep(0, n), sigma = Sigma_new, log = TRUE) +
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



