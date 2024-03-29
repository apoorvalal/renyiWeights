% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebweights.R
\name{entrBal}
\alias{entrBal}
\title{Compute entropy balancing weights that minimize KL divergence between uniform and solution weights. These are Renyi weights for α = 1}
\usage{
entrBal(
  target_moments,
  X,
  max.iterations = 200L,
  base.weights = NULL,
  constraint.tolerance = 1,
  print.level = 0
)
}
\arguments{
\item{target_moments}{[vector, matrix] target moments to match X columns (1 x k)}

\item{X}{[data.frame, matrix] table of covariates for source sample (n x k)}

\item{max.iterations}{[200] stopping rule}

\item{base.weights}{= [NULL] n-vector of baseline weights}

\item{constraint.tolerance}{[1] value for constraint threshold}

\item{print.level}{[0, 1, 2, 3] 0 is silent, 1 reports success, 2 and 3 are noisy (for debugging)}
}
\value{
list containing a n-vector of weights
}
\description{
Compute entropy balancing weights that minimize KL divergence between uniform and solution weights. These are Renyi weights for α = 1
}
