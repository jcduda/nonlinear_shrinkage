
reduce_res <- function(res, method, true_fn, samples=NULL){

  if(length(res) == 1 & all(is.na(res)))
    return(list(RMSE = NA, y_hat = NA, upper = NA, lower = NA))

  #-----------------------------------------------------------------------------

  if(method %in% c("nlfs", "bspline", "parametric_bspline", "pspline")) {
    all_y <- res$draws$all_y_hat[, samples]
    doses <- res$data$x
    resp <- res$data$y
  }

  if(method == "parametric"){
    all_y <- predict(res)$Y
    doses <- res$data[, 1]
    resp <- res$data[, 2]
  }

  #-----------------------------------------------------------------------------

  y_hat <- rowMeans(all_y)
  y_true <- true_fn(doses)

  upper <- apply(all_y, 1, function(x) quantile(x, 0.975))
  lower <- apply(all_y, 1, function(x) quantile(x, 0.025))

  RMSE <- sqrt(mean( (y_true - y_hat)^2 ))

  return(list(RMSE = RMSE,
              y_hat = y_hat,
              upper = upper,
              lower = lower,
              x = doses,
              y = resp))
}
