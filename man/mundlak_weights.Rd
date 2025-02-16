% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mundlak.R
\name{mundlak_weights}
\alias{mundlak_weights}
\title{Re-weight control sub-groups to treated sub-group means}
\usage{
mundlak_weights(
  ind_covs,
  interact_covs,
  trt,
  Z,
  lambda = 0,
  lowlim = 0,
  uplim = 1,
  scale_sample_size = TRUE,
  exact_global = TRUE,
  verbose = TRUE,
  eps_abs = 1e-05,
  eps_rel = 1e-05,
  ...
)
}
\arguments{
\item{ind_covs}{n x d1 matrix of covariates for individual units}

\item{interact_covs}{n x d2 matrix of interactions between individual and cluster covariates}

\item{trt}{Vector of treatment assignments}

\item{Z}{Vector of group indicators with J levels}

\item{lambda}{Regularization hyper parameter, default 0}

\item{lowlim}{Lower limit on weights, default 0}

\item{uplim}{Upper limit on weights, default 1}

\item{scale_sample_size}{Whether to scale the dispersion penalty by the sample size of each group, default TRUE}

\item{exact_global}{Whether to enforce exact balance for overall population on individual covariates}

\item{verbose}{Whether to show messages, default T}

\item{eps_abs}{Absolute error tolerance for solver}

\item{eps_rel}{Relative error tolerance for solver}

\item{...}{Extra arguments for osqp solver}
}
\value{
\itemize{
         \item{weights }{Estimated weights as a length n vector}
         \item{imbalance }{Imbalance in covariates as a d X J matrix}
         \item{global_imbalance}{Overall imbalance in covariates, as a length d vector }}
}
\description{
Re-weight control sub-groups to treated sub-group means
}
