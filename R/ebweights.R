# %% ####################################################
#' Compute entropy balancing weights that minimize KL divergence between uniform and solution weights. These are Renyi weights for α = 1
#' @param X [data.frame, matrix] table of covariates for source sample (n x k)
#' @param target_moments [vector, matrix] target moments to match X columns (1 x k)
#' @param  max.iterations [200] stopping rule
#' @param base.weights = [NULL] n-vector of baseline weights
#' @param constraint.tolerance [1] value for constraint threshold
#' @param print.level [0, 1, 2, 3] 0 is silent, 1 reports success, 2 and 3 are noisy (for debugging)
#' @return list containing a n-vector of weights
#' @export
entrBal = function(X, target_moments,
                   max.iterations = 200L,
                   base.weights = NULL,
                   constraint.tolerance = 1, print.level = 0) {
  X = as.matrix(X);
  nsource = nrow(X)
  if (is.null(base.weights)) base.weight = rep(1, nsource)
  # n0 X C
  X0 = cbind(1, X)
  # rank check
  if (qr(X0)$rank != ncol(X0))
    stop("collinearity in covariate matrix for controls (remove collinear covariates)")
  # intercept plus target moments
  # C-vector
  X1m = c(1, target_moments)
  # default coefficients - C-vector of 0s
  coefs = rep(0, ncol(X0))
  ## run algo
  eb.out = eb_solve_dual(
    X1m                  = X1m,
    X0                   = X0,
    coefs                = coefs,
    base.weight          = base.weight,
    max.iterations       = max.iterations,
    constraint.tolerance = constraint.tolerance,
    print.level          = print.level
  )
  if (eb.out$converged == TRUE & print.level > 0)
    cat("Converged within tolerance \n")
  z = eb.out$Weights.ebal
  return(z)
}

#' Entropy balancing by solving primal
#' Mainly for didactic purposes; don't use.
#' @param X1m K-vector of target means
#' @param X0 NXK matrix of covariates
#' @param solv Solver for CVXR ()
#' @import CVXR
#' @export
eb_solve_primal = function(X1m, X0, solv = "MOSEK") {
  wt = Variable(nrow(X0)) # weights are of length n_0
  objective = Maximize(sum(entr(wt))) # entropy objective fn
  constraints = list(
    wt >= 0, sum(wt) == 1, # proper weights
    t(X0) %*% wt == X1m # balance
  )
  prob = Problem(objective, constraints)
  result = solve(prob, solver = solv)
  wt_hat = result$getValue(wt)
  return(wt_hat)
}


#' Entropy balancing by solving dual
#' @param X1m K-vector of target means
#' @param X0 NXK matrix of covariates
#' @param base.weight = [NULL] n-vector of baseline weights
#' @param max.iterations [200] stopping rule
#' @param constraint.tolerance [1] value for constraint threshold
#' @param print.level [0, 1, 2, 3] 0 is silent, 1 reports success, 2 and 3 are noisy (for debugging)
#' @param sparsify [T/F] (in progress) run Newton-Raphson with sparse matrix classes from Matrix package
#' @import Matrix
#' @export
eb_solve_dual = function(X1m, X0,
                         coefs = NULL, base.weight = NULL,
                         max.iterations = 200L,
                         constraint.tolerance = 1,
                         print.level = 0,
                         sparisfy = FALSE) {
  if (is.null(coefs)) coefs = rep(0, ncol(X0))
  if (is.null(base.weight)) base.weight = rep(1, nrow(X0))

  # !TODO make this work with sparse matrices
  # if (sparsify) {
  #   X0          = as(X0, "dgCMatrix")
  #   base.weight = as(base.weight, "sparseVector")
  #   X1m         = as(X1m, "sparseVector")
  # }

  converged = FALSE
  for (iter in 1:max.iterations) {
    # Z is a R-vector of solution coefficients (coefs)
    #########################################################
    # unnormalised weights
    weights.temp = c(exp(X0 %*% coefs))
    weights.ebal = weights.temp * base.weight
    ### gradient construction
    X0.agg = c(weights.ebal %*% X0)
    # ∇_Z = M - CW
    gradient = X0.agg - X1m
    # minimum reached
    if (max(abs(gradient)) < constraint.tolerance) {
      converged = TRUE
      break
    }
    # noisy
    if (print.level >= 2)
      cat(
        "Iteration", iter, "maximum deviation is =",
        format(max(abs(gradient)), digits = 4), "\n"
      )
    hessian = t(X0) %*% (weights.ebal * X0)
    Coefs = coefs
    newton = solve(hessian, gradient)
    coefs = coefs - newton

    # step length - optimal step length section
    loss.new = line.searcher(
      Base.weight = base.weight, X0 = X0,
      X1m = X1m, coefs = coefs, Newton = newton, ss = 1
    )
    loss.old = line.searcher(
      Base.weight = base.weight, X0 = X0,
      X1m = X1m, coefs = Coefs, Newton = newton, ss = 0
    )

    if (is.na(loss.new) | is.na(loss.old)) {
      stop(
        "Optimization ran into problems. Loss function values are",
        "\n New:", loss.new,
        "\n Old:", loss.old,
        "\n Constraint tolerance is", constraint.tolerance,
        "\n Increase constraint tolerance or trim moment conditions problem"
      )
    }

    if (print.level >= 3) cat("new loss", loss.new, "old loss=", loss.old, "\n")
    if (loss.old <= loss.new) {
      ss.out = optimize(line.searcher,
        lower = .00001, upper = 1, maximum = FALSE,
        Base.weight = base.weight, X0 = X0, X1m = X1m,
        coefs = Coefs, Newton = newton
      )
      if (print.level >= 3)
        cat("LS Step Length is ", ss.out$minimum, "\n")
      if (print.level >= 3)
        cat("Loss is", ss.out$objective, "\n")
      coefs = Coefs - ss.out$minimum * solve(hessian, gradient)
    }
  }
  # step out of loop
  if (print.level >= 1 && converged) cat("Converged within tolerance \n")

  return(
    list(
      maxdiff = max(abs(gradient)),
      coefs = coefs,
      Weights.ebal = weights.ebal,
      converged = converged
    )
  )
}

# %% internal for step length
line.searcher = function(Base.weight, X0, X1m, coefs, Newton, ss) {
  weights.temp = c(exp(X0 %*% (coefs - (ss * Newton))))
  weights.temp = weights.temp * Base.weight
  X0.agg = c(weights.temp %*% X0)
  maxdiff = max(abs(X0.agg - X1m))
  return(maxdiff)
}
