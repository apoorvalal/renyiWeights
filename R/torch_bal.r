# %% ####################################################
#' Balancing weights using dual formulation implemented in torch (autodiff via torch)
#' @param X0 matrix of donor covariates
#' @param X1 target moments
#' @param lr learning rate for optimizer
#' @param num_iterations number of iterations to run the optimizer
#' @param objective loss function (either "entr" or "l2")
#' @param noi boolean - print last value of gradient
#' @param debug print lots of intermediate data
#' @return vector of weights
#' @export
#' @import torch

bal_lbfgs = function(
    X0, X1, lr = 0.01, num_iterations = 50,
    objective = c("entr", "l2"), noi = FALSE, debug = FALSE) {
  obj = match.arg(objective)
  # list inputs
  nc = nrow(X0)
  inp = list(x0 = torch_tensor(X0), x1 = torch_tensor(X1))
  ######################################################################
  # loss fns
  ebal_loss = \(lambda) {
    inner_sum = torch_exp(-1 * torch_matmul(inp$x0, lambda))
    loss_value = torch_log(torch_sum(inner_sum)) + lambda$dot(inp$x1)
    loss_value
  }
  l2bal_loss = \(lambda) {
    xb = torch_matmul(inp$x0, lambda)
    inner_sum = (-1 * xb^2) / 4 + xb / nc
    loss_value = -1 * torch_sum(inner_sum) + lambda$dot(inp$x1)
    loss_value
  }
  loss = switch(obj,
    l2 = l2bal_loss,
    entr = ebal_loss
  )
  ######################################################################
  # coefficients starting values
  lam = torch_rand(ncol(X0), requires_grad = TRUE)
  # optimizer
  optimizer = optim_lbfgs(lam, lr = lr, line_search_fn = "strong_wolfe")
  # anonymous function that lbfgs uses
  calc_loss = \() {
    optimizer$zero_grad()
    value = loss(lam)
    if (debug) cat("Value is", as.numeric(value), "\n")
    value$backward()
    value
  }
  # iterate
  for (i in 1:num_iterations) {
    optimizer$step(calc_loss)
    if (debug) cat("iter ", as.numeric(i), "\n")
    if (debug) cat("gradient:", as.numeric(lam$grad), "\n")
  }
  if (noi) cat("Final gradient:", as.numeric(lam$grad), "\n")
  # take final value of lam
  λ = as.numeric(lam)
  # weights from lagrangian
  wts = switch(obj,
    l2 = (-1 * as.matrix(inp$x0) %*% λ) / 2 + 1 / nc,
    entr = exp(-1 * as.matrix(inp$x0) %*% λ)
  )
  # normalised
  wts = as.numeric(wts / sum(wts))
}

# %%
#' ebal implementation with autodiff via torch
#' @param X0 donor units matrix
#' @param X1 target moments
#' @return vector of weights
#' @export
ebal_torch = function(X0, X1) {
  inp = list(x0 = torch_tensor(X0), x1 = torch_tensor(X1))
  # loss fn returns tensor
  ebal_loss_torch = \(lambda) {
    inner_sum = torch_exp(-1 * torch_matmul(inp$x0, lambda))
    loss_value = torch_log(torch_sum(inner_sum)) + torch_matmul(lambda, inp$x1)
  }
  # gradient - autograd
  loss_grad = \(lambda) {
    lambda_ad = torch_tensor(lambda, requires_grad = TRUE)
    loss = ebal_loss_torch(lambda_ad)
    grad = autograd_grad(loss, lambda_ad)[[1]]
    as.numeric(grad)
  }
  # call optim - fn needs to have numeric output
  λ = optim(
    fn = \(x) as.numeric(ebal_loss_torch(x)),
    gr = loss_grad, par = rep(1, ncol(X0)),
    method = "BFGS"
  )$par
  # extract weights from lagrangian
  wts = exp(-1 * as.matrix(inp$x0) %*% λ)
  # normalised
  wts / sum(wts)
}
