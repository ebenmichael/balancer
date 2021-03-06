% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{balancer_boot}
\alias{balancer_boot}
\title{Get point estimate and bootstraped standard errors}
\usage{
balancer_boot(
  nboot,
  X,
  trt,
  Z = NULL,
  type = c("att", "subgrp", "missing", "hte"),
  link = c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
  regularizer = c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc"),
  hyperparam,
  Q = NULL,
  kernel = NULL,
  kern_param = 1,
  normalized = TRUE,
  sample_trt = FALSE,
  opts = list()
)
}
\arguments{
\item{nboot}{Number of bootstrap samples}

\item{X}{n x d matrix of covariates}

\item{trt}{Vector of treatment status indicators}

\item{Z}{Vector of subgroup indicators or observed indicators}

\item{type}{Find balancing weights for ATT, subgroup ATTs,
ATT with missing outcomes, and heterogeneouts effects}

\item{link}{Link function for weights}

\item{regularizer}{Dual of balance criterion}

\item{hyperparam}{Regularization hyperparameter}

\item{Q}{m x m matrix to tie together ridge penalty, default: NULL,
if TRUE, use covariance of treated groups}

\item{kernel}{What kernel to use, default NULL}

\item{kern_param}{Hyperparameter for kernel}

\item{normalized}{Whether to normalize the weights}

\item{sample_trt}{Whether to resample treated units}

\item{opts}{Optimization options
\itemize{
  \item{MAX_ITERS }{Maximum number of iterations to run}
  \item{EPS }{Error rolerance}
  \item{alpha }{Elastic net parameter}}}
}
\value{
\itemize{
         \item{theta }{Estimated dual propensity score parameters}
         \item{weights }{Estimated primal weights}
         \item{imbalance }{Imbalance in covariates}}
}
\description{
Get point estimate and bootstraped standard errors
}
