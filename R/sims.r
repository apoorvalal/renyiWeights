# %%
#' DGP from Frolich (2007) / Hainmueller (2012)
#' @param error_design string (1, 2, or 3) for low, medium, or high overlap.
#' @param outcome_design string (1, 2, or 3) for linear, complex, or quadratic outcome model.
#' @param n sample size
#' @export
hainmueller = function(error_design, outcome_design, n = 1e4) {
  ###################################################
  # generate X
  ###################################################
  # covariance matrix
  Σ = matr(
    2, 1, -1 |
      1, 1, -.5 |
      -1, -.5, 1
  )
  X13 = MASS::mvrnorm(n = n, rep(0, 3), Σ)
  X4 = runif(n = n, -3, 3)
  X5 = rchisq(n = n, 1)
  X6 = rbinom(n = n, 1, 0.5)
  X = cbind(X13, X4, X5, X6)
  colnames(X) = paste0("X_", 1:6)
  ###################################################
  # generate W
  ###################################################
  # generate error term
  gen_ε = \(d){
    ε = switch(d,
      "1" = rnorm(n = n, 0, sqrt(30)),
      "2" = 0.5 + (rchisq(n, df = 5) - 5) * sqrt(67.6 / 10), # scale var by 10 because vχ2 = 2k
      "3" = rnorm(n, 0, sqrt(100))
    )
  }
  ε = gen_ε(error_design)
  βw = c(1, 2, -2, -1, -.5, 1)
  # generate treatment
  w = 1 * (X %*% βw + ε > 0)
  ###################################################
  # generate Y
  ###################################################
  gen_y = \(d, X){
    y = switch(d,
      "1" = X[, 1] + X[, 2] + X[, 3] - X[, 4] + X[, 5] + X[, 6],
      "2" = X[, 1] + X[, 2] + 0.2 * X[, 3] * X[, 4] + sqrt(X[, 5]),
      "3" = (X[, 1] + X[, 2] + X[, 5])^2
    )
  }
  y = gen_y(outcome_design, X) + rnorm(n)
  # return
  return(
    list(y = y, w = w, X = X)
  )
}
