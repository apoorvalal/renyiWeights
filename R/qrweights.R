#' Balancing weights to minimise quadratic divergence between uniform and solution weights. These are Renyi weights for α = 2
#' @param X [data.frame, matrix] table of covariates for source sample (n x k)
#' @param target_moments [vector, matrix] target moments to match X columns (1 x k)
#' @param base weights n_0 vector of base weights
#' @param solv solver for CVXR
#' @return list containing a n-vector of weights
#' @import CVXR
#' @export
qr_solve_primal = function(X1m, X0, base_weights = NULL, solv = "MOSEK") {
  wt = Variable(nrow(X0)) # weights are of length n_0
  if (is.null(base_weights)) {
    objective = Minimize(sum(wt^2)) # quadratic weights
  } else {
    objective = Minimize(sum(wt^2 / base_weights)) # quadratic weights
  }
  constraints = list(
    sum(wt) == 1, # now allowed to be negative, just need to sum to 1
    t(X0) %*% wt == X1m # balance
  )
  prob = Problem(objective, constraints)
  result = solve(prob, solver = 'MOSEK')
  wt_hat = result$getValue(wt)
}

# %%
#' Compute OLS with arbitrary vector of weights (possibly negative)
#' @param y Response vector
#' @param X A numeric data matrix
#' @param w vector of weights (same length as y)
#' @return Regression vector beta of length ncol(X).
#' @export
#' @import car
OLSw = function(y, X, w) {
  XtWX = car::wcrossprod(X, w = w)
  XtWy = car::wcrossprod(X, y, w = w)
  solve(XtWX) %*% XtWy
}
# %%
