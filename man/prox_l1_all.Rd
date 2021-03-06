% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{prox_l1_all}
\alias{prox_l1_all}
\title{L1 Prox for multilevel model, separate for global/local params + intercepts}
\usage{
prox_l1_all(x, lam, opts)
}
\arguments{
\item{x}{Input matrix (two sets of parameters x = U + V)}

\item{lam}{Prox scaling factor}

\item{opts}{List of options (opts["lam"] holds the other scaling}
}
\value{
soft thresholded parameters with different soft thresholds
}
\description{
L1 Prox for multilevel model, separate for global/local params + intercepts
}
