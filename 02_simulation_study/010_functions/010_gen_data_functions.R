library(DoseFinding)
library(dplyr)


gen_data <- function(true_fn, sigma_sq, n){
  x <- sort(runif(n))
  y <- true_fn(x) + rnorm(n, sd = sqrt(sigma_sq))
  d <- data.frame(x = x, y = y)
  res <- list(x = x,
              y = y,
              true_fn = true_fn,
              n = n,
              sigma_sq = sigma_sq)
}


get_true_fn <- function(truth){

  # For hill and hilldown
  if(truth %in% c("hill", "hilldown")){

    all_non_lin_theta <- c(ED50 = 0.3, h = 6)

    true_fn <- function(x) sigEmax(x, e0 = 0, eMax = 1,
                                   ed50 = all_non_lin_theta["ED50"],
                                   h = all_non_lin_theta["h"])


    if(truth == "hilldown"){
      hill_part <- true_fn
      down_part <- function(x) ifelse(x < 0.6, 0, -1.5*(x-0.6)^2)
      true_fn <- function(x) hill_part(x) + down_part(x)
    }

    attributes(true_fn)$theta <- all_non_lin_theta
  }

  if(truth == "power"){
    true_fn <- function(x) (0 + 1*x^0.5)
    attributes(true_fn)$theta = c(d = 0.5)
  }

  return(true_fn)
}

#get_n_knots <- function(n) round(sqrt(n)) # use fixed n_knots = 15

get_prior <- function(model){
  if(model == "hill"){
      theta <- c(ED50_mean = 0.5, ED50_var = 0.05, h_mean = 3, h_var = 3)
  }

  if(model == "power"){
      theta <- c(d_mean = 0.5, d_var = 0.25, NA, NA)
  }

  if(model == "hill_power"){
    theta <- c(ED50_mean = 0.5, ED50_var = 0.05, h_mean = 3, h_var = 3,
               d_mean = 0.5, d_var = 0.25)
  }

  return(theta)
}

