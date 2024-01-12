nlfs <- function(n_draw = 3000, n_inner_knots = 10, d,
                 nlf_deriv, theta_start = NULL, true_fn = NULL,
                 assume = "hill",
                 shrinkage = "own_slice", # or: "half_cauchy"
                 # Add smoothing in B:
                 add_smooth=FALSE,
                 # hyperparams for omega
                 omega_a = 0.5,
                 omega_b = NULL, # defined based on number of knots
                 # hyper params for sigma^2~IG(sq_a, sq_b)
                 sigma_sq_a = 0.001,
                 sigma_sq_b = 0.001,
                 # hyper params for intercept ~ N(mu_interc, sq_interc)
                 interc_mu = 0,
                 interc_sq = 1,
                 update_theta = TRUE,
                 # hyper params for ED50 (theta 3) ~ N(mu_ed50, sq_ed50)
                 theta3_mu = 0.5,
                 theta3_sq = 0.05,
                 # hyper params for h (theta 4)
                 theta4_mean = 3,
                 theta4_var = 3,
                 # if assume=hill_power: last theta is for power: d
                 theta5_mu = NA,
                 theta5_sq = NA,
                 # extended for, e.g. testosterone data that is not in (0,1)
                 sigma_sq_half_cauchy = FALSE,
                 lower = NULL,
                 upper = NULL,
                 boundary_knots = NULL
){



  n <- d$n
  # Derivative functions -------------------------------------------------------
  if(assume == "hill"){
    s <- 2 # dimension of non-lin parameters
    if(is.null(theta_start)) theta_start <- c(0.5, 5)
    nlf_deriv <- function(x, theta){
      stopifnot(length(theta)==2)
      sigEmaxGrad(x, 1, theta[1], theta[2])
    }
  }
  if(assume == "power"){
    s <- 1
    if(is.null(theta_start)) theta_start <- c(1.2)
    nlf_deriv <- function(x, theta){
      stopifnot(length(theta)== 1)
      # from DoseFinding::sigEmaxGrad
      lg2 <- function(x) {
        l <- x
        l[x == 0] <- 0
        l[x != 0] <- log(x[x != 0])
        l
      }
      cbind(1, b=x^theta[1], d=lg2(x)*x^theta[1]) # keep intercept!
    }
  }

  if(assume == "hill_power"){
    s <- 3
    if(is.null(theta_start)) theta_start <- c(0.5, 5, 1.2)
    nlf_deriv <- function(x, theta){
      stopifnot(length(theta)== 3)
      # first 4 columns from hill_deriv:

      # from DoseFinding::sigEmaxGrad
      lg2 <- function(x) {
        l <- x
        l[x == 0] <- 0
        l[x != 0] <- log(x[x != 0])
        l
      }
      cbind(sigEmaxGrad(x, 1, theta[1], theta[2]), # b=x^theta[3],
            d=lg2(x)*x^theta[3]) # keep intercept!
    }
  }
  # prepare B ------------------------------------------------------------------
  if(is.null(lower)) lower <- min(d$x)
  if(is.null(upper)) upper <- max(d$x)
  if(is.null(boundary_knots)) boundary_knots <- c(0,1)
  prepare_B <- prep_B(d$x, n_inner_knots, lower = lower, upper = upper, boundary_knots = boundary_knots)
  B <- prepare_B$B

  k <- ncol(B)
  R <- diag(k) # could be abandoned. Was once used to incorporate smoothing.


  # set omega_b as in Shin
  if(is.null(omega_b)) omega_b <- exp(-k*log(n) / 2)

  # prepare containers ---------------------------------------------------------

  # per draw: one column
  all_beta <- matrix(0, ncol = n_draw, nrow = k)
  all_omega <- rep(0, n_draw)
  all_tau_sq <- rep(0, n_draw) # always consider both parametrizations
  all_xi <- rep(0, n_draw) # for shrinkage = "half_cauchy"
  all_sigma_sq <- rep(0, n_draw)
  all_xi_sigma <- rep(0, n_draw) # for testosterone data (sigma_sq_half_cauchy == TRUE)
  all_interc <- rep(0, n_draw)
  all_theta <- matrix(0, ncol = n_draw, nrow = s)
  all_y_hat <- matrix(0, ncol = n_draw, nrow = n)

  # Initialize -----------------------------------------------------------------

  all_interc[1] <- mean(d$y[1:round(0.1*n)])
  all_omega[1] <- 0.5
  all_tau_sq[1] <- 1 / all_omega[1] - 1
  all_xi[1] <- 0.5
  all_xi_sigma[1] <- 0.5
  all_beta[, 1] <- coef(lm(d$y-all_interc[1] ~ B - 1))
  all_y_hat[, 1] <- all_interc[1] + B %*% all_beta[, 1]
  all_sigma_sq[1] <- var(d$y - all_y_hat[, 1])
  all_theta[, 1] <- theta_start

  # Updating ---------------------------------------------------------------------

  # calculate only once
  sigma_sq_a_new <- (n + k) / 2 + sigma_sq_a
  if(grepl("hill", assume)){
    h_log_var <- log(theta4_var / theta4_mean^2 + 1)
    h_log_mean <- log(theta4_mean) - h_log_var/2
  }

  set.seed(1234)
  for(i in 2:n_draw){
    if(i %% 50 == 0) print(i)
    #print(i)
    F_mat <- nlf_deriv(d$x, all_theta[, i-1])
    # Use no more:
    # P <- F_mat %*% solve(t(F_mat) %*% F_mat) %*% t(F_mat)
    # Use qr decomposition instead
    L <- t(suppressWarnings(chol(crossprod(F_mat), pivot = TRUE)))
    d0 <- attr(L, "rank")
    Qt <- forwardsolve(L, t(F_mat[, attr(L, "pivot")]), d0)
    P <- crossprod(Qt)

    BIPB <- t(R) %*% t(B) %*% (diag(n) - P) %*% B %*% R

    # beta update
    all_beta[, i] <- beta_update(y = d$y, # no group-index for y-tilde
                                 Sigma_prec = BIPB,
                                 tau_sq = all_tau_sq[i-1],
                                 sigma_sq = all_sigma_sq[i-1],
                                 alpha = all_interc[i-1],
                                 X_mat = B)

    # intercept update
    all_interc[i] <- interc_update(y = d$y, X_mat = B, n = n, alpha_mu = interc_mu,
                                   alpha_v = interc_sq, sigma_sq = all_sigma_sq[i-1],
                                   beta = all_beta[, i])

    # sigma_sq update (+ extension for testosterone: sigma^2~half_cuachy for overdispersion)
    all_sigma_sq[i] <- sigma_sq_update(sigma_sq_a_new = sigma_sq_a_new,
                                       sigma_sq_b = ifelse(sigma_sq_half_cauchy == FALSE, sigma_sq_b, 1/all_xi_sigma[i-1]),
                                       n = n, k = k, tau_sq = all_tau_sq[i-1], beta = all_beta[, i],
                                       Sigma_prec = BIPB, alpha = all_interc[i], X_mat = B, y = d$y)

    if(sigma_sq_half_cauchy == TRUE){
      all_xi_sigma[i] <- rinvgamma(1, 1, 1 + 1 / all_sigma_sq[i])
    }



    # tau^2 / omega update
    # all_tau_sq[i] <- 1e-16
    if(shrinkage == "own_slice"){
      all_tau_sq[i] <- tau_sq_update(beta = all_beta[, i], old_tau_sq = all_tau_sq[i-1],
                                     omega_a = omega_a, omega_b = omega_b,
                                     k = k, Sigma_prec = BIPB,
                                     sigma_sq = all_sigma_sq[i])
    }

    if(shrinkage == "half_cauchy"){
      all_tau_sq[i] <- tau_sq_update_hc(beta = all_beta[, i], Sigma_prec = BIPB,
                                        xi = all_xi[i-1], sigma_sq = all_sigma_sq[i-1])
      all_xi[i] <- rinvgamma(1, 1, 1 + 1 / all_tau_sq[i])
    }

    all_omega[i] <- 1 / (1 +  all_tau_sq[i])


    # update non-linear parameters!

    #set.seed(4)
    if(update_theta == TRUE){

      XPX_curr_inv <- chol2inv(chol(BIPB))
      Sigma_curr <- all_sigma_sq[i] * all_tau_sq[i] * B %*% XPX_curr_inv %*% t(B) + all_sigma_sq[i] * diag(n)

      if(grepl("hill", assume)){

        # Theta 3 (ED50 in hill)
        # used for all non-lin params:
        all_theta[1, i] <- ED50_update(y = d$y - all_interc[i],  tau_sq = all_tau_sq[i],
                                       sigma_sq = all_sigma_sq[i],  ED50_curr = all_theta[1, i-1],
                                       h = all_theta[2, i-1],
                                       d_par = ifelse(assume == "hill_power", all_theta[3, i-1], NA),
                                       x_grid = d$x,
                                       X = B, Sigma_curr = Sigma_curr,
                                       nlf_deriv = nlf_deriv,
                                       ED50_mu = theta3_mu,
                                       ED50_v = theta3_sq,
                                       ED50_v_prop = ifelse(i < 100, theta3_sq, var(all_theta[1, (i-99):i])))

        # Theta 4 (h or steepness in hill)

        # if(i >=100){
        #   h_var <- var(all_theta[2, (i-99):i])
        #   h_mean <- mean(all_theta[2, (i-99):i])
        #   adap_h_log_var <- log(h_var / h_mean^2 + 1)
        # }

        all_theta[2, i] <- h_update(y = d$y - all_interc[i], sigma_sq = all_sigma_sq[i],
                                    ED50 = all_theta[1, i],
                                    d_par = ifelse(assume == "hill_power", all_theta[3, i-1], NA),
                                    h_curr = all_theta[2, i-1],
                                    x_grid = d$x, X = B, Sigma_curr = Sigma_curr,
                                    nlf_deriv = nlf_deriv,
                                    h_mu = h_log_mean,
                                    h_v = h_log_var,
                                    h_v_prop = h_log_var,#ifelse(i<100, h_log_var, adap_h_log_var),
                                    tau_sq = all_tau_sq[i])
      }


      if(grepl("power", assume)){
        # Similar to ED50 update: non-linear parameter with assumed normal prior
        # in power-model, the only non-linear parameter is the d

        if(assume == "hill_power") d_pos <- 3
        if(assume == "power") d_pos <- 1

        all_theta[d_pos, i] <- d_update(y = d$y - all_interc[i],  tau_sq = all_tau_sq[i],
                                    sigma_sq = all_sigma_sq[i],
                                    d_curr = ifelse(grepl("hill", assume), all_theta[d_pos, i-1], all_theta[1, i-1]),
                                    ED50 = ifelse(grepl("hill", assume), all_theta[1, i-1], NA),
                                    h = ifelse(grepl("hill", assume), all_theta[2, i-1], NA),
                                    x_grid = d$x,
                                    nlf_deriv = nlf_deriv,
                                    X = B, Sigma_curr = Sigma_curr,
                                    d_mu = ifelse(grepl("hill", assume), theta5_mu, theta3_mu),
                                    d_v = ifelse(grepl("hill", assume), theta5_sq, theta3_sq),
                                    d_v_prop = ifelse(i<100, ifelse(grepl("hill", assume), theta5_sq, theta3_sq), var(all_theta[d_pos, (i-99):i])))

      }
    }

    # The actual fit:
    all_y_hat[, i] <- all_interc[i] + B %*% all_beta[, i]

  }

  return(res = list(data = d,
                    B = B,
                    R = R,
                    inner_knots = prepare_B$inner_knots,
                    n_final_inner_knots = prepare_B$final_n_inner_knots,

                    inputs = list(n = n,
                                  n_draw = n_draw,
                                  n_inner_knots = n_inner_knots,
                                  k = k,
                                  add_smooth = add_smooth,
                                  nlf_deriv = nlf_deriv,
                                  shrinkage = shrinkage,
                                  update_theta = update_theta,
                                  theta_start = theta_start,
                                  true_fn = true_fn,
                                  assume = assume,
                                  omega_a = omega_a,
                                  omega_b = omega_b,
                                  sigma_sq_a = sigma_sq_a,
                                  sigma_sq_b = sigma_sq_b,
                                  sigma_sq_half_cauchy = sigma_sq_half_cauchy,
                                  interc_mu = interc_mu,
                                  interc_sq = interc_sq,
                                  theta3_mu = theta3_mu,
                                  theta3_sq = theta3_sq,
                                  theta4_mean = theta4_mean,
                                  theta4_var = theta4_var,
                                  theta5_mu = theta5_mu,
                                  theta5_sq = theta5_sq),
                    draws = list(all_beta = all_beta,
                                 all_omega = all_omega,
                                 all_tau_sq = all_tau_sq,
                                 all_sigma_sq = all_sigma_sq,
                                 all_interc = all_interc,
                                 all_theta = all_theta,
                                 all_y_hat = all_y_hat)))


}
