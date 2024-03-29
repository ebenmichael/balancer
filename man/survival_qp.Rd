% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_qp.R
\name{survival_qp}
\alias{survival_qp}
\title{Re-weight control sub-groups to treated sub-group means}
\usage{
survival_qp(
  B_X,
  trt,
  times,
  events,
  t,
  lambda = 0,
  lowlim = 0,
  uplim = NULL,
  verbose = TRUE,
  eps_abs = 1e-05,
  eps_rel = 1e-05,
  ...
)
}
\arguments{
\item{B_X}{n x k basis matrix of covariates}

\item{trt}{Vector of treatment assignments}

\item{times}{Vector of event/censoring times (see events for if event or censoring time)}

\item{events}{Vector of boolean censoring indicators (whether individual was censored)}

\item{t}{Time}

\item{lambda}{Regularization hyperparameter, default 0}

\item{lowlim}{Lower limit on weights, default 1}

\item{uplim}{Upper limit on weights, default NULL}

\item{verbose}{Whether to show messages, default T}

\item{eps_abs}{Absolute error tolerance for solver}

\item{eps_rel}{Relative error tolerance for solver}

\item{...}{Extra arguments for osqp solver}
}
\value{
\itemize{
         \item{weights}{Estimated weights as a length n vector}
         \item{imbalance0}{Imbalance in covariates for control in vector}
         \item(imbalance1){Imbalance in covariates for treated in vector}}
}
\description{
Re-weight control sub-groups to treated sub-group means
}
