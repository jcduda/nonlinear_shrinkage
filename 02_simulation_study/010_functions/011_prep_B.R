prep_B <- function(x, n_inner_knots, lower = min(x), upper = max(x), boundary_knots = c(0, 1)){
  knots_ok <- 0
  final_n_inner_knots <- n_inner_knots

  while(knots_ok == 0){
    inner_knots <- seq(lower + 0.5*1/final_n_inner_knots, upper - 0.5*1/final_n_inner_knots,
                       length.out = n_inner_knots)
    B <- bs(x = x,
            degree = 3,
            knots = inner_knots,
            Boundary.knots = boundary_knots,
            intercept = F)

    if(any(colSums(B)==0)){
      final_n_inner_knots <- final_n_inner_knots - 1
      warning(paste0("At least one column in B was zero. Redo B with ",
                     final_n_inner_knots, " inner knots instead of original ",
                     n_inner_knots, " inner knots."))

      inner_knots <- seq(lower + 0.5*1/final_n_inner_knots, upper - 0.5*1/final_n_inner_knots,
                         length.out = final_n_inner_knots)
      B <- bs(x = x,
              degree = 3,
              knots = inner_knots,
              Boundary.knots = boundary_knots,
              intercept = F)

    } else {
      knots_ok <- 1
    }
  }

  return(list(B = B,
         n_inner_knots = n_inner_knots,
         final_n_inner_knots = final_n_inner_knots,
         inner_knots = inner_knots))
}
