# Bayesian splines (no smoothing), simply: B-splines

# Bayesian P-spline implementation

bspline <- function(n_draw = 3000, n_inner_knots = 10, d = d,
                     sigma_sq_a = 0.001,
                     sigma_sq_b = 0.001,
                    lambda_sq_a = 0.001,
                    lambda_sq_b = 0.001,
                     # hyper params for intercept ~ N(mu_interc, sq_interc)
                     interc_mu = 0,
                     interc_sq = 1
){

  ################################################################################
  n <- d$n


  # prepare the B-spline matrix
  prepare_B <- prep_B(d$x, n_inner_knots)
  B <- prepare_B$B
  k <- ncol(B)

  # beta ~ N(0, sigma^2 I)
  # Precision = K^(-1) = solve(t(R) %*% R)

  prec <- diag(k)

  # prepare containers -----------------------------------------------------------

  # per draw: one column
  all_beta <- matrix(0, ncol = n_draw, nrow = k)
  all_sigma_sq <- rep(0, n_draw)
  all_lambda_sq <- rep(0, n_draw) # IG-scaling prior for beta covariance
  all_interc <- rep(0, n_draw)
  all_y_hat <- matrix(0, ncol = n_draw, nrow = n)

  # Initialize

  all_interc[1] <- mean(d$y[1:round(0.1*n)])
  all_beta[, 1] <- coef(lm(d$y-all_interc[1] ~ B - 1))
  all_lambda_sq[1] <- var(all_beta[, 1])
  all_y_hat[, 1] <- all_interc[1] + B %*% all_beta[, 1]
  all_sigma_sq[1] <- var(d$y - all_y_hat[, 1])


  #plot(d$x, B%*%all_beta[, 1], type = "l")

  # Updating ---------------------------------------------------------------------

  # calculate only once
  sigma_sq_a_new <- (n + k) / 2 + sigma_sq_a

  lambda_sq_a_new <- k/2 + lambda_sq_a

  for(i in 2:n_draw){
    if(i %% 500 == 0) print(i)
    # beta update
    all_beta[, i] <- beta_update(y = d$y, # no group-index for y-tilde
                                 Sigma_prec = prec,
                                 tau_sq = all_lambda_sq[i-1],
                                 sigma_sq = all_sigma_sq[i-1],
                                 alpha = all_interc[i-1],
                                 X_mat = B)

    # intercept update
    all_interc[i] <- interc_update(y = d$y, X_mat = B, n = n, alpha_mu = interc_mu,
                                   alpha_v = interc_sq, sigma_sq = all_sigma_sq[i-1],
                                   beta = all_beta[, i])

    # sigma_sq update
    all_sigma_sq[i] <- sigma_sq_update(sigma_sq_a_new = sigma_sq_a_new, sigma_sq_b = sigma_sq_b,
                                       n = n, k = k, tau_sq = all_lambda_sq[i-1], beta = all_beta[, i],
                                       Sigma_prec = prec, alpha = all_interc[i], X_mat = B, y = d$y)

    # lambda_sq update
    all_lambda_sq[i] <- rinvgamma(1, shape = lambda_sq_a_new,
                                  rate = 1/(2*all_sigma_sq[i]) * crossprod(all_beta[, i]) + lambda_sq_b)


    # The actual fit:
    all_y_hat[, i] <- all_interc[i] + B %*% all_beta[, i]
  }

  return(res = list(data = d,
                    B = B,
                    inner_knots = prepare_B$inner_knots,
                    n_final_inner_knots = prepare_B$final_n_inner_knots,
                    inputs = list(n = n,
                                  n_draw = n_draw,
                                  n_inner_knots = n_inner_knots,
                                  k = k,
                                  sigma_sq_a = sigma_sq_a,
                                  sigma_sq_b = sigma_sq_b,
                                  lambda_sq_a = lambda_sq_a,
                                  lambda_sq_b = lambda_sq_b,
                                  interc_mu = interc_mu,
                                  interc_sq = interc_sq),
                    draws = list(all_beta = all_beta,
                                 all_sigma_sq = all_sigma_sq,
                                 all_lambda_sq = all_lambda_sq,
                                 all_interc = all_interc,
                                 all_y_hat = all_y_hat)))


}

