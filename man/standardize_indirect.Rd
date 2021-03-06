% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{standardize_indirect}
\alias{standardize_indirect}
\title{Re-weight populations to group targets}
\usage{
standardize_indirect(
  X,
  Z,
  lambda = 0,
  lowlim = 0,
  uplim = 1,
  scale_sample_size = F,
  verbose = TRUE,
  n_cores = 1,
  eps_abs = 1e-05,
  eps_rel = 1e-05,
  ...
)
}
\arguments{
\item{X}{n x d matrix of covariates}

\item{Z}{Vector of group indicators with J levels}

\item{lambda}{Regularization hyper parameter, default 0}

\item{lowlim}{Lower limit on weights, default 0}

\item{uplim}{Upper limit on weights, default 1}

\item{scale_sample_size}{Whether to scale the dispersion penalty by the sample size of each group, default T}

\item{verbose}{Whether to show messages, default T}

\item{n_cores}{Number of cores to find weights in parallel}

\item{eps_abs}{Absolute error tolerance for solver}

\item{eps_rel}{Relative error tolerance for solver}

\item{...}{Extra arguments for osqp solver}
}
\value{
\itemize{
         \item{weights }{Estimated weights as an n x J matrix}
         \item{imbalance }{Imbalance in covariates as a d X J matrix}
         }
}
\description{
Re-weight populations to group targets
}
