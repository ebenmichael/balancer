% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_treatment.R
\name{create_constraints_treatment}
\alias{create_constraints_treatment}
\title{Create the constraints for QP: l <= Ax <= u}
\usage{
create_constraints_treatment(
  X0s,
  Xtaus,
  target,
  Z,
  S_factor,
  pro_trt_split,
  pro_ctr_split,
  lowlim,
  uplim,
  verbose
)
}
\arguments{
\item{X0s}{List of J n x d0 matrices of covariates split by group}

\item{Xtaus}{List of J n x dtau matrices of covariates split by group}

\item{target}{Vector of population means to re-weight to}

\item{Z}{Vector of group indicators}

\item{S_factor}{Vector of site indicators}

\item{pro_trt_split}{List of J vectors of treatment propensity multipliers}

\item{pro_ctr_split}{List of J vectors of control propensity multipliers}

\item{lowlim}{Lower limit on weights}

\item{uplim}{Upper limit on weights}

\item{verbose}{Boolean indicating whether to print progress}
}
\value{
A, l, and u
}
\description{
Create the constraints for QP: l <= Ax <= u
}
