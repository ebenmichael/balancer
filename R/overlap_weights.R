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
  return(list(weights = wts, imbalance = imbal))
}

create_constraints_maxsize <- function(X, trt, lowlim, uplim, verbose, exact_balance) {

  n <- nrow(X)
  d <- ncol(X)
  n1 <- sum(trt)
  n0 <- sum(1 - trt)

  # sum to number of treated/control units
  A1 <- rbind(trt, (1 - trt))
  A1 <- cbind(A1, Matrix::Matrix(0, 2, d))
  l1 <- c(n1, n0)
  u1 <- c(n1, n0)


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
  m <- length(unique(clusters))
  q <- numeric(n + d + m)
  P <- create_P_matrix_cluster_overlap(trt, d, m, lambda , icc, FALSE)


  constraints <- create_constraints_maxsize_cluster(X, trt, clusters, lowlim, uplim, verbose, FALSE)

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


create_P_matrix_cluster_overlap <- function(trt, d, m, lambda, icc, exact_balance) {

  n <- length(trt)
  # first include a diagonal element for the auxiliary covariates X %*% gamma
  P1 <- Matrix::bdiag(Matrix::Matrix(0, n, n),
                      Matrix::Diagonal(d),
                      Matrix::Matrix(0, m, m))
  # Add iid variance term
  I0 <- Matrix::bdiag(Matrix::Diagonal(n) * (trt / sum(trt)^2 + (1 - trt) / sum(1 - trt)^2),
                      Matrix::Diagonal(d, 0),
                      Matrix::Matrix(0, m, m))
  P2 <- lambda * (1 - icc) * I0

  # add correlation within cluster term
  P3 <- lambda * icc *  Matrix::bdiag(Matrix::Diagonal(n, 0),
                                     Matrix::Diagonal(d, 0),
                                     Matrix::Diagonal(m))

  P <- P1 + P2 + P3
  return(P)
}




create_constraints_maxsize_cluster <- function(X, trt, clusters, lowlim, uplim, verbose, exact_balance) {

  n <- nrow(X)
  d <- ncol(X)
  n1 <- sum(trt)
  n0 <- sum(1 - trt)
  cluster_mat <- Matrix::sparse.model.matrix(~ as.factor(clusters) - 1) *
                    (trt / sum(trt) + (1 - trt) / sum(1 - trt))
  m <- ncol(cluster_mat)

  # sum to number of treated/control units
  A1 <- rbind(trt, (1 - trt))
  A1 <- cbind(A1, Matrix::Matrix(0, 2, d))
  A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow = nrow(A1), ncol = m))
  l1 <- c(n1, n0)
  u1 <- c(n1, n0)


  # upper and lower bounds
  A2 <-  Matrix::Diagonal(n)
  A2 <- cbind(A2, Matrix::Matrix(0, n, d))
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = m))
  l2 <- lowlim * (trt * n1  + (1 - trt) * n0)
  u2 <- uplim * (trt * n1  + (1 - trt) * n0)

  # auxiliary variable
  Xtrt <- (trt / sum(trt) - (1 - trt) / sum(1 - trt)) * X
  A3 <- cbind(t(Xtrt), -diag(d))
  A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = m))
  l3 <- numeric(d)
  u3 <- numeric(d)

  # create intermediate variable that is the (normalized) sum of weights in a cluster
  A4 <- Matrix::cbind2(Matrix::t(cluster_mat), Matrix::Matrix(0, nrow = m, ncol = d))
  A4 <- Matrix::cbind2(A4, -Matrix::Diagonal(m))
  l4 <- rep(0, m)
  u4 <- rep(0, m)

  # if exact_balance is true, then constrain weights to exactly balance treatment and control
  if(exact_balance) {
    A5 <- cbind(matrix(0,ncol = n, nrow = d), diag(d))
    A5 <- Matrix::cbind2(A5, Matrix::Matrix(0, nrow = nrow(A5), ncol = m))
    l5 <- numeric(d)
    u5 <- numeric(d)

    A <- rbind(A1, A2, A3, A4, A5)
    l <- c(l1, l2, l3, l4, l5)
    u <- c(u1, u2, u3, u4, u5)
  } else {
    A <- rbind(A1, A2, A3, A4)
    l <- c(l1, l2, l3, l4)
    u <- c(u1, u2, u3, u4)
  }

  return(list(A = A, l = l, u = u))
}



#' Compute point estimate and standard error with clustered max subset weights
#' @param y Vector of outcomes
#' @param wts Vector of weight
#' @param trt Vector of treatment assignments
#' @param clusters Vector of cluster assignments
#' @param m1hat Vector of model predictions of E[Y(1) | covariates]
#' @param m0hat Vector of model predictions of E[Y(0) | covariates]
#' 
#' @return Data.Frame with the point estimate and standard error
#' @export
compute_cluster_maxsubset_se <- function(y, wts, trt, clusters, m1hat, m0hat) {

  cluster_mat <- Matrix::sparse.model.matrix(~ as.factor(clusters) - 1)
  n1 <- sum(trt)
  n0 <- sum(1 - trt)


  mu1 <- sum(m1hat * wts) / sum(wts) +  sum((y - m1hat)[trt == 1] * wts[trt == 1]) / sum(wts[trt == 1])
  mu0 <- sum(m0hat * wts) / sum(wts) +  sum((y - m0hat)[trt == 0] * wts[trt == 0]) / sum(wts[trt == 0])

  mhat <- m1hat * trt + m0hat * (1 - trt)
  resids <- y - mhat


  wtd_resids1 <- Matrix::t(cluster_mat[trt == 1, ]) %*% (wts * resids)[trt == 1]
  wtd_resids0 <- Matrix::t(cluster_mat[trt == 0, ]) %*% (wts * resids)[trt == 0]
  wtd_resids <- Matrix::t(cluster_mat) %*% (wts * resids)
  m1hat_clus <- Matrix::t(cluster_mat[trt == 1, ]) %*% (wts[trt == 1] * (m1hat[trt == 1] - mu1))
  m0hat_clus <- Matrix::t(cluster_mat[trt == 0, ]) %*% (wts[trt == 0] * (m0hat[trt == 0] - mu0))
  tauhat_clus <- Matrix::t(cluster_mat) %*% (wts * ((m1hat - m0hat) - (mu1 - mu0)))
  se1 <- sqrt(sum(wtd_resids1^2)/sum(wts[trt == 1])^2 + sum(m1hat_clus^2) / sum(wts[trt == 1])^2)
  se0 <- sqrt(sum(wtd_resids0^2)/sum(wts[trt == 0])^2 +  sum(m0hat_clus^2) / sum(wts[trt == 0])^2)
  se_tau <- sqrt(sum(wtd_resids^2)/sum(wts)^2 +  sum(tauhat_clus^2) / sum(wts)^2)


  return(data.frame(Estimand = c("Average Outcome for Treated", "Average Counterfactual Outcome for Treated", "ATT"),
                    Estimate = c(mu1, mu0, mu1 - mu0), SE = c(se1, se0, se_tau)))
}