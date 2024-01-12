parametric_bspline <- function(d = d,
                               n_draw = 6000, n_inner_knots = 10,
                               sigma_sq_a = 0.001,
                               sigma_sq_b = 0.001,
                               # We shrink B-spline betas using shin shrinkage
                               omega_a = 0.5,
                               omega_b = NULL, # defined based on number of knots
                               # hyper params for intercept ~ N(mu_interc, sq_interc)
                               interc_mu = 0,
                               interc_sq = 1,
                               assume = "hill",
                               shrinkage = "half_cauchy", # own_slice or half_cauchy
                               #theta_start,
                               true_fn = NULL,
                               # hyper params for ED50 (theta 3) ~ N(mu_ed50, sq_ed50)
                               # Or for d in case of power model
                               theta3_mu = 0.5,
                               theta3_sq = 0.05,
                               # hyper params for h (theta 4)
                               theta4_mean = 3,
                               theta4_var = 3,
                               # In param_bspline also: prior for Emax:
                               emax_mu = 1.5,
                               emax_sq = 2){

  n <- d$n

  # Count non-linear parameters
  if(assume == "hill") {
    s <- 2
  }
  if(assume == "power"){
    s <- 1
  }

  if(!(assume %in% c("hill", "power")))
    stop("Function only for assume = hill or assume = power")
  # Our model: y = hill(theta) + B-spline + error,
  #            or power(theta)


  # prepare the B-spline matrix
  prepare_B <- prep_B(d$x, n_inner_knots)
  B <- prepare_B$B

  k <- ncol(B)

  # beta ~ N(0, sigma^2 tau^2 diag(all_lambda_sq))
  prec <- diag(k)

  # set omega_b as in Shin
  if(is.null(omega_b)) omega_b <- exp(-k*log(n) / 2)

  # prepare containers ---------------------------------------------------------

  # per draw: one column
  all_beta <- matrix(0, ncol = n_draw, nrow = k)
  all_omega <- rep(0, n_draw)
  all_tau_sq <- rep(0, n_draw) # always consider both parametrizations
  all_xi <- rep(0, n_draw) # for shrinkage = "half_cauchy"

  # Don't forget: vanilla horse shoe also has local shrinkage!
  all_lambda_sq <- matrix(0, ncol = n_draw, nrow = k)
  all_nu <- matrix(0, ncol = n_draw, nrow = k)

  all_sigma_sq <- rep(0, n_draw)
  all_interc <- rep(0, n_draw)
  # Only non-linear
  all_theta <- matrix(0, ncol = n_draw, nrow = s)
  # Emax parameter in hill or scalingparameter b for power model
  # note: in this function, param_bspline, we also model all linear parameters
  all_emax <- rep(0, n_draw)
  all_y_hat <- matrix(0, ncol = n_draw, nrow = n)

  # calculate only once
  sigma_sq_a_new <- (n + k) / 2 + sigma_sq_a
  if(assume == "hill"){
    h_log_var <- log(theta4_var / theta4_mean^2 + 1)
    h_log_mean <- log(theta4_mean) - h_log_var/2
  }

  # Initialize -----------------------------------------------------------------


  all_omega[1] <- 0.5
  all_tau_sq[1] <- 1 / all_omega[1] - 1
  all_xi[1] <- 0.5

  all_nu[, 1] <- rinvgamma(k, 1/2, 1)
  all_lambda_sq[, 1] <- rinvgamma(k, 1/2, 1/all_nu[, 1])

  if(!(assume %in% c("hill", "power")))
    stop("param_bspline so far only implemented for assume=hill or assume=power")
  # y = h(theta) + b-spline + error, model intercept separately

  if(assume %in% "hill"){
    curr_coef <- coef(fitMod(dose = d$x, resp = d$y - all_interc[1],
                             model = "sigEmax", type = "normal"))
    all_interc[1] <- curr_coef[1]
    all_emax[1] <- curr_coef[2]
    all_theta[, 1] <- curr_coef[3:4]

    curr_param_part <- sigEmax(d$x, e0 = 0,
                               eMax = all_emax[1],
                               ed50 = all_theta[1, 1],
                               h = all_theta[2, 1])
  }

  if(assume %in% "power"){
    # Just to get starting values
    power_fit <- single_continuous_fit(d$x,
                                       d$y,
                                       model_type = "power",
                                       BMR_TYPE="sd", BMR = 0.1,
                                       #burnin = burnin,
                                       samples = 10000
    )
    all_interc[1] <- power_fit$parameters[1]
    all_emax[1] <- power_fit$parameters[2]
    all_theta[, 1] <- power_fit$parameters[3]

    power_model <- function(x, a, b, d) a + b*x^d

    curr_param_part <- power_model(d$x, a=0,
                                   b = all_emax[1],
                                   d = all_theta[1, 1])

  }


  all_beta[, 1] <- coef(lm(d$y- all_interc[1] - curr_param_part  ~ B - 1))
  all_y_hat[, 1] <- all_interc[1] + curr_param_part + B %*% all_beta[, 1]
  all_sigma_sq[1] <- var(d$y - all_y_hat[, 1])

  set.seed(1234)
  for(i in 2:n_draw){
    if(i %% 500 == 0) print(i)

    prec <- diag(1 / all_lambda_sq[, i-1])
    # beta_update
    all_beta[, i] <- beta_update(y = d$y, # no group-index for y-tilde
                                 Sigma_prec = prec,
                                 tau_sq = all_tau_sq[i-1],
                                 sigma_sq = all_sigma_sq[i-1],
                                 alpha = all_interc[i-1] + curr_param_part,
                                 X_mat = B)

    # intercept update
    all_interc[i] <- interc_update(y = d$y - curr_param_part, # is this not allowed!?
                                   X_mat = B, n = n, alpha_mu = interc_mu,
                                   alpha_v = interc_sq, sigma_sq = all_sigma_sq[i-1],
                                   beta = all_beta[, i])

    # sigma_sq update
    all_sigma_sq[i] <- sigma_sq_update(sigma_sq_a_new = sigma_sq_a_new, sigma_sq_b = sigma_sq_b,
                                       n = n, k = k, tau_sq = all_tau_sq[i-1], beta = all_beta[, i],
                                       Sigma_prec = prec, alpha = all_interc[i], X_mat = B,
                                       y = d$y - curr_param_part)
    # tau_sq update
    # (always use own shin implementation with omega and a and b)

    if(shrinkage == "own_slice"){
      stop("For parametric_b_spline only (vanilla) horseshoe with half cauchy")
      all_tau_sq[i] <- tau_sq_update(beta = all_beta[, i], old_tau_sq = all_tau_sq[i-1],
                                     omega_a = omega_a, omega_b = omega_b,
                                     k = k, Sigma_prec = prec,
                                     sigma_sq = all_sigma_sq[i])
    }


    if(shrinkage == "half_cauchy"){
      all_tau_sq[i] <- tau_sq_update_hc(beta = all_beta[, i], Sigma_prec = prec,
                                        xi = all_xi[i-1], sigma_sq = all_sigma_sq[i])
      all_xi[i] <- rinvgamma(1, 1, 1 + 1 / all_tau_sq[i])

      all_lambda_sq[, i] <- rinvgamma(k, 1,
                                     1/all_nu[, i-1] + (all_beta[, i])^2 / (2*all_tau_sq[i]*all_sigma_sq[i]))

      all_nu[, i] <- rinvgamma(k, 1, 1 + 1/all_lambda_sq[, i])
    }

    all_omega[i] <- 1 / (1 +  all_tau_sq[i])

    # Update parameters from the parametric part

    # Linear Emax parameter from hill:
    if(assume == "hill"){
      h_theta0 <- sigEmax(d$x, 0, 1, all_theta[1, i-1], all_theta[2, i-1])
    }
    if(assume == "power"){
      h_theta0 <- power_model(d$x, 0, 1, all_theta[1, i-1])
    }
    all_emax[i] <- emax_update(y_tilde = d$y - all_interc[i] - B %*% all_beta[, i],
                               h_theta0 = h_theta0,
                               sigma_sq = all_sigma_sq[i],
                               emax_sq = emax_sq, emax_mu = emax_mu)

    # Update non-linear parameters -------------------------------------------

    if(assume == "hill"){ #---------------------------------------------------

      # ED50
      # For some reason: super low acceptance rate...
      all_theta[1, i] <- update_ED50_param_bspline(x = d$x,
                                                   ED50_curr = all_theta[1, i-1],
                                                   h = all_theta[2, i-1],
                                                   ED50_mu = theta3_mu,
                                                   ED50_v_prop = ifelse(i<100, theta3_sq, var(all_theta[1, (i-99):i])),
                                                   ED50_v = theta3_sq,
                                                   y_tilde = (d$y - all_interc[i] - B %*% all_beta[, i]) / all_emax[i],
                                                   sigma_sq_tilde =  all_sigma_sq[i] / (all_emax[i])^2)

      # h
      # if(i >=100){
      #   h_var <- var(all_theta[2, (i-99):(i-1)])
      #   h_mean <- mean(all_theta[2, (i-99):(i-1)])
      #   adap_h_log_var <- log(h_var / h_mean^2 + 1)
      # }

      all_theta[2, i] <- update_h_param_bspline(x = d$x,
                                                h_curr = all_theta[2, i-1],
                                                ED50 = all_theta[1, i],
                                                h_v = h_log_var,
                                                h_v_prop = h_log_var, #ifelse(i<100, h_log_var, adap_h_log_var),

                                                h_mu = h_log_mean,
                                                y_tilde = (d$y - all_interc[i] - B %*% all_beta[, i]) / all_emax[i],
                                                sigma_sq_tilde = all_sigma_sq[i] / (all_emax[i])^2)

      curr_param_part <-  sigEmax(d$x, e0 = 0, eMax = all_emax[i], ed50 = all_theta[1, i], h = all_theta[2, i])
    }
    if(assume == "power"){ #---------------------------------------------------

      all_theta[1, i] <- update_d_param_bspline(x = d$x,
                                                d_curr = all_theta[1, i-1],
                                                d_mu = theta3_mu,
                                                d_v_prop = ifelse(i < 100, theta3_sq, var(all_theta[1, (i-99):i])),
                                                d_v = theta3_sq,
                                                y_tilde = (d$y - all_interc[i] - B %*% all_beta[, i]) / all_emax[i],
                                                sigma_sq_tilde = all_sigma_sq[i] / (all_emax[i])^2)

      curr_param_part <-  power_model(d$x, a = 0, b = all_emax[i], d = all_theta[1, i])

    }


    # The actual fit
    # The actual fit:
    all_y_hat[, i] <- all_interc[i] + B %*% all_beta[, i] +
      # parametric part:
      curr_param_part
  }

  return(res = list(data = d,
                    B = B,
                    inner_knots = prepare_B$inner_knots,
                    n_final_inner_knots = prepare_B$final_n_inner_knots,

                    inputs = list(n = n,
                                  n_draw = n_draw,
                                  n_inner_knots = n_inner_knots,
                                  k = k,
                                  # nlf_deriv = nlf_deriv,
                                  shrinkage = shrinkage,
                                  # update_theta = update_theta,
                                  # theta_start = theta_start,
                                  true_fn = true_fn,
                                  assume = assume,
                                  omega_a = omega_a,
                                  omega_b = omega_b,
                                  # sigma_sq_half_cauchy = sigma_sq_half_cauchy,
                                  sigma_sq_a = sigma_sq_a,
                                  sigma_sq_b = sigma_sq_b,
                                  interc_mu = interc_mu,
                                  interc_sq = interc_sq,
                                  theta3_mu = theta3_mu,
                                  theta3_sq = theta3_sq,
                                  theta4_mean = theta4_mean,
                                  theta4_var = theta4_var,
                                  emax_mu = emax_mu,
                                  emax_sq = emax_sq), # continue here
                    draws = list(all_beta = all_beta,
                                 all_omega = all_omega,
                                 all_tau_sq = all_tau_sq,
                                 all_sigma_sq = all_sigma_sq,
                                 all_lambda_sq = all_lambda_sq,
                                 all_nu = all_nu,
                                 all_interc = all_interc,
                                 all_theta = all_theta,
                                 all_emax = all_emax,
                                 all_y_hat = all_y_hat)))

}