#'  Balancing weights using pogs
#' @param X0 matrix of donor covariates
#' @param X1 target moments
#' @param loss one of "entr" or "l2"
#' @param tol tolerance for balance constraint violation
#' @return vector of weights
#' @export
#' @import pogs
bal_pogs = function(X0, X1, loss = c("entr", "l2"), tol = 0){
  loss = match.arg(loss)
  # loss switcher
  lossfn = switch(loss,
    entr = kNegEntr(),
    l2   = kSquare()
  )
  # pogs form of balancing problem
  # Minimize ||gamma||_2^2 , ST
  # ||X0'gamma - X1||_infty <= tol and Σ gamma_i = 1.
  # Equivalently, in notation recognized by POGS, we can write
  # Minimize ||gamma||_2^2 +
  # subject to I[1:2p] <= 0, J[2p+1, 2p+2] = 1, where
  #     (X0'  -X1 -tol)   (gamma)
  # I = (-X0' X1  -tol) * (ONE  )
  #     (1'  0  0      )   (    )
  #     (0   1  0      )
  g = list(h = lossfn)
  f = list(
    h = c(kIndLe0(2 * ncol(X0)), kIndEq0(2)),
    b = c(rep(0, 2 * ncol(X0)), 1, 1)
  )
  A = rbind(
    cbind(t(X0), -X1, -tol),
    cbind(-t(X0), X1, -tol),
    c(rep(1, nrow(X0)), 0, 0),
    c(rep(0, nrow(X0)), 1, 0)
  )
  f$h = c(f$h, kIndLe0(nrow(X0)))
  f$b = c(f$b, rep(0, nrow(X0)))
  A = rbind(A, cbind(diag(-1, nrow(X0)), 0, 0))
  # solve
  pogs.solution = pogs(A, f, g, params = list(rel_tol = 1e-4, abs_tol = 1e-5))
  gamma = pogs.solution$x[1:nrow(X0)]
  ω_pogs = gamma / sum(gamma)
}
