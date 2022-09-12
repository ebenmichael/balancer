#' Find maximum effective sample size balanced set
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment assignments
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1 * number of treated/control units
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
maxsubset_weights <- function(X, trt, lambda = 0, lowlim = 0, uplim = 1,
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

create_constraints_maxsize <- function(X, trt, lowlim, uplim, verbose, exact_balance) {

  n <- nrow(X)
  d <- ncol(X)
  n1 <- sum(trt)
  n0 <- sum(1 - trt)

  # sum to number of treated/control units
  A1 <- rbind(trt / n1, (1 - trt) / n0)
  A1 <- cbind(A1, Matrix::Matrix(0, 2, d))
  l1 <- c(1, 1)
  u1 <- c(1, 1)


  # upper and lower bounds
  A2 <-  Matrix::Diagonal(n)
  A2 <- cbind(A2, Matrix::Matrix(0, n, d))
  l2 <- lowlim * (trt * n1  + (1 - trt) * n0)
  u2 <- uplim * (trt * n1  + (1 - trt) * n0)

  # auxiliary variable
  Xtrt <- (trt / sum(trt) - (1 - trt) / sum(1 - trt)) * X
  A3 <- cbind(t(Xtrt), -diag(d))
  l3 <- numeric(d)
  u3 <- numeric(d)

  # if exact_balance is true, then constrain weights to exactly balance treatment and control
  if(exact_balance) {
    A4 <- cbind(matrix(0,ncol = n, nrow = d), diag(d))
    l4 <- numeric(d)
    u4 <- numeric(d)
    print(dim(A4))

    A <- rbind(A1, A2, A3, A4)
    l <- c(l1, l2, l3, l4)
    u <- c(u1, u2, u3, u4)
  } else {
    A <- rbind(A1, A2, A3)
    l <- c(l1, l2, l3)
    u <- c(u1, u2, u3)
  }

  return(list(A = A, l = l, u = u))
}


# #' Find maximum effective sample size balanced set for clustered observational studies
# #' @param ind_covs n x d1 matrix of covariates for individual units
# #' @param clus_covs n x d2 matrix of covariates for clusters
# #' @param trt Vector of treatment assignments
# #' @param lambda Regularization hyper parameter, default 0
# #' @param lowlim Lower limit on weights, default 0
# #' @param uplim Upper limit on weights, default 1 1 * number of treated/control units
# #' @param verbose Whether to show messages, default T
# #' @param eps_abs Absolute error tolerance for solver
# #' @param eps_rel Relative error tolerance for solver
# #' @param ... Extra arguments for osqp solver
# #'
# #' @return \itemize{
# #'          \item{weights }{Estimated weights as a length n vector}
# #'          \item{imbalance }{Imbalance in covariates as a d X J matrix}
# #'          \item{global_imbalance}{Overall imbalance in covariates, as a length d vector }}
# #' @export
# maxsubset_weights_cluster <- function(ind_covs, clus_covs, trt, lambda = 0,
#                                       lowlim = 0, uplim = 1, verbose = TRUE,
#                                       eps_abs = 1e-5, eps_rel = 1e-5, ...) {

#   X <- cbind(ind_covs, clus_covs)

#   out <- maxsubset_weights(X, trt, lambda, lowlim, uplim, verbose, eps_abs, eps_rel, ...)

#   out$ind_imbalance <- out$imbalance[1:ncol(ind_covs)]
#   out$clus_imbalance = out$imbalance[(ncol(ind_covs) + 1):(ncol(X))]
#   return(out[-2])
# }


#' Find maximum effective sample size balanced set
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param clus_covs n x d2 matrix of covariates for clusters
#' @param trt Vector of treatment assignments
#' @param clusters n vector of cluster assignments
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param icc Intraclass correlation coefficient for regularization
#' @param uplim Upper limit on weights, default 1 * number of treated/control units
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
maxsubset_weights_cluster <- function(ind_covs, clus_covs, trt, clusters,
                                      lambda = 0, icc = 0, lowlim = 0, uplim = 1,
                                      verbose = TRUE, eps_abs = 1e-5, eps_rel = 1e-5,
                                      ...) {
  
  # convert X to a matrix
  X <- cbind(ind_covs, clus_covs)
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  q <- numeric(n + d)
  P <- create_P_matrix_cluster_overlap(trt, d, clusters, lambda , icc, FALSE)


  constraints <- create_constraints_maxsize(X, trt, lowlim, uplim, verbose, FALSE)

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
  ind_imbalance <- imbal[1:ncol(ind_covs)]
  clus_imbalance = imbal[(ncol(ind_covs) + 1):(ncol(X))]
  return(list(weights = wts,
              ind_imbalance = ind_imbalance,
              clus_imbalance = clus_imbalance))
}


create_P_matrix_cluster_overlap <- function(trt, d, clusters, lambda, icc, exact_balance) {

  n <- length(trt)
  # first include a diagonal element for the auxiliary covariates X %*% gamma
  P1 <- Matrix::bdiag(Matrix::Matrix(0, n, n),
                        Matrix::Diagonal(d))
  # Add iid variance term
  I0 <- Matrix::bdiag(Matrix::Diagonal(n) * (trt / sum(trt)^2 + (1 - trt) / sum(1 - trt)^2),
                      Matrix::Diagonal(d, 0))
  P2 <- lambda * (1 - icc) * I0

  # add correlation within cluster term
  cluster_mat <- Matrix::sparse.model.matrix(~ as.factor(clusters) - 1) *
                    (trt / sum(trt) + (1 - trt) / sum(1 - trt))
  P3 <- Matrix::bdiag(lambda * icc * cluster_mat %*% Matrix::t(cluster_mat),
                      Matrix::Diagonal(d, 0))

  P <- P1 + P2 + P3
  return(P)
}


