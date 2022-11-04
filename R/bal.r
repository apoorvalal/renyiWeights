# %%
#' Bal: Approximate balancing weights
#'
#' @description Solve for a set of balancing weights subject to a convex loss.
#'
#'
#' @param Xt          Target means (K vector)
#' @param X           Data matrix (N × K)
#' @param objective   Loss function for weights. One of entropy, L1, L2.
#' @param tol         Tolerance for imbalance. 0 by default, which corresponds with exact balance.
#' @param solv        Solver name to pass to CVXR solve
#' @return n-vector of weights
#' @export

bal = function(Xt, X,
               objective = c("entropy", "l1", "l2"),
               tol = 0,
               solv = "MOSEK") {
  # loss
  objective = match.arg(objective)
  # preprocessing
  X = as.matrix(X); n = nrow(X); k = ncol(X)
  # optim
  ω = Variable(n, nonneg = TRUE)
  δ = Parameter(1, nonneg = TRUE)
  obj = switch(objective,
    l1      = p_norm((ω - 1 / n), p = 1),
    l2      = p_norm((ω - 1 / n), p = 2),
    entropy = -sum(entr(ω))
  )
  if (is.null(obj)) error("Invalid objective function")
  objective = Minimize(obj) # entropy objective fn
  constraints = list(
    ω >= 0,
    sum(ω) == 1, # proper weights
    abs(t(X) %*% ω) - Xt <= δ # balance within tolerance
  )
  value(δ) = tol
  prob = Problem(objective, constraints)
  result = solve(prob, solver = solv)
  ω_hat = result$getValue(ω)
  return(as.numeric(ω_hat))
}
