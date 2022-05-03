#' Find maximum effective sample size balanced set
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment assignments
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param verbose Whether to show messages, default T
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#'
#' @return \itemize{
#'          \item{weights }{Estimated weights as a length n vector}
#'          \item{imbalance }{Imbalance in covariates as a d X J matrix}
#'          \item{global_imbalance}{Overall imbalance in covariates, as a length d vector }}
#' @export
maxsubset_weights <- function(X, trt, lambda = 0, lowlim = 0, uplim = Inf,
                              verbose = TRUE, eps_abs = 1e-5, eps_rel = 1e-5,
                              ...) {
  # convert X to a matrix
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  q <- c(-rep(1, n) * lambda / n, numeric(d))
  P <- Matrix::bdiag(Matrix::Diagonal(n, lambda / n), Matrix::Diagonal(d))


  constraints <- create_constraints_maxsize(X, trt, lowlim, uplim, verbose)

  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                              eps_rel = eps_rel,
                              eps_abs = eps_abs),
                      list(...)))

  solution <- osqp::solve_osqp(P, q, constraints$A,
                                  constraints$l, constraints$u,
                                  pars = settings)

  wts <- solution$x[1:n]
  imbal <- colSums((trt / sum(trt) - (1  - trt) / sum(1 - trt)) * wts * X)
  return(list(weights = wts, imbalance = imbal))
}

create_constraints_maxsize <- function(X, trt, lowlim, uplim, verbose) {

  n <- nrow(X)
  d <- ncol(X)

  # sum to number of treated/control units
  A1 <- rbind(trt / sum(trt), (1 - trt) / sum(1 - trt))
  A1 <- cbind(A1, Matrix::Matrix(0, 2, d))
  l1 <- c(1, 1)
  u1 <- c(1, 1)


  # upper and lower bounds
  A2 <-  Matrix::Diagonal(n)
  A2 <- cbind(A2, Matrix::Matrix(0, n, d))
  l2 <- rep(lowlim, n)
  u2 <- rep(uplim, n)

  # auxiliary variable
  Xtrt <- (trt / sum(trt) - (1 - trt) / sum(1 - trt)) * X
  A3 <- cbind(t(Xtrt), -diag(d))
  l3 <- numeric(d)
  u3 <- numeric(d)

  A <- rbind(A1, A2, A3)
  l <- c(l1, l2, l3)
  u <- c(u1, u2, u3)

  return(list(A = A, l = l, u = u))
}