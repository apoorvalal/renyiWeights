% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebweights.R
\name{eb_solve_dual}
\alias{eb_solve_dual}
\title{Entropy balancing by solving dual}
\usage{
eb_solve_dual(
  X1m,
  X0,
  coefs = NULL,
  base.weight = NULL,
  max.iterations = 200L,
  constraint.tolerance = 1,
  print.level = 0,
  sparisfy = FALSE
)
}
\arguments{
\item{X1m}{K-vector of target means}

\item{X0}{NXK matrix of covariates}

\item{base.weight}{= [NULL] n-vector of baseline weights}

\item{max.iterations}{[200] stopping rule}

\item{constraint.tolerance}{[1] value for constraint threshold}

\item{print.level}{[0, 1, 2, 3] 0 is silent, 1 reports success, 2 and 3 are noisy (for debugging)}

\item{sparsify}{[T/F] (in progress) run Newton-Raphson with sparse matrix classes from Matrix package}
}
\description{
Entropy balancing by solving dual
}
