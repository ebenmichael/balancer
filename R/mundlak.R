################################################################################
## Balancing weights for clustered data using cluster-level sufficient stats
## as in Mundlak regression
################################################################################

#' Re-weight control sub-groups to treated sub-group means
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param interact_covs n x d2 matrix of interactions between individual and cluster covariates
#' @param trt Vector of treatment assignments
#' @param Z Vector of group indicators with J levels
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default TRUE
#' @param exact_global Whether to enforce exact balance for overall population on individual covariates
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
mundlak_weights <- function(ind_covs, interact_covs, trt, Z, lambda = 0, lowlim = 0, uplim = 1,
                            scale_sample_size = TRUE, exact_global = TRUE,
                            verbose = TRUE,
                            eps_abs = 1e-5, eps_rel = 1e-5, ...) {

  # convert ind_covs and interact_covs to a matrix
  ind_covs <- as.matrix(ind_covs)
  interact_covs <- as.matrix(interact_covs)
  d_ind <- ncol(ind_covs)
  d_interact <- ncol(interact_covs)
  aux_dim <- d_ind + d_interact


  # split data and treatment by factor
  Z_factor <- as.factor(Z)
  Xz <- split.data.frame(ind_covs, Z_factor)
  trtz <- split(trt, Z)

  check_data_multi(ind_covs, trt, Z, Xz, lambda, lowlim, uplim)


  unique_Z <- levels(Z_factor)
  J <- length(unique_Z)
  n <- nrow(ind_covs)




  idxs <- split(1:nrow(ind_covs), Z_factor)



  # Setup the components of the QP and solve
  if(verbose) message("Creating linear term vector...")
  q <- create_q_vector_mundlak(ind_covs, interact_covs, trt)

  if(verbose) message("Creating quadratic term matrix...")
  P <- create_P_matrix_multi(n, aux_dim)

  I0 <- create_I0_matrix_multi(Xz, scale_sample_size, n, aux_dim)
  P <- P + lambda * I0

  if(verbose) message("Creating constraint matrix...")
  constraints <- create_constraints_mundlak(ind_covs, interact_covs, trt, Z, lowlim,
                                            uplim, exact_global, verbose)

  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                              eps_rel = eps_rel,
                              eps_abs = eps_abs),
                      list(...)))

  solution <- osqp::solve_osqp(P, q, constraints$A,
                                  constraints$l, constraints$u,
                                  pars = settings)

  weights <- solution$x[1:n]
  # compute imbalance matrix

  # compute individual covariates imbalance
  global_imbal <- t(ind_covs) %*% weights / sum(weights) - colMeans(ind_covs[trt == 1,, drop = F])
  interact_imbal <- t(interact_covs) %*% weights / sum(weights) - colMeans(interact_covs[trt == 1,, drop = F])
  return(list(weights = weights,
              interaction_imbalance = interact_imbal,
              global_imbalance = global_imbal))

}

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param interact_covs n x d2 matrix of interactions between individual and cluster covariates
#'
#' @return q vector
create_q_vector_mundlak <- function(ind_covs, interact_covs, trt) {
  aux_dim <- ncol(ind_covs) + ncol(interact_covs)
  n <- nrow(ind_covs)
  q <- - c(colSums(ind_covs[trt == 1,, drop = F]),
           colSums(interact_covs[trt == 1,, drop = F])
          )

  q <- Matrix::sparseVector(q, (n + 1):(n + aux_dim),
                            n + aux_dim)
  return(q)
}


#' Create the constraints for QP: l <= Ax <= u
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param interact_covs n x d2 matrix of interactions between individual and cluster covariates
#' @param trt Vector of treatment assignments
#' @param Z Vector of group indicators with J levels
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param exact_global Boolean whether to include an exact global constraint
#' @param verbose Boolean whether to display progress
#'
#' @return A, l, and u
create_constraints_mundlak <- function(ind_covs, interact_covs, trt, Z, lowlim,
                                     uplim, exact_global, verbose) {

  n <- nrow(ind_covs)

  # dimension of auxiliary weights
  aux_dim <- ncol(ind_covs) + ncol(interact_covs)


  if(verbose) message("\tx Sum to one constraint")
  # sum-to-n1j constraint for each group
  A1_small <- Matrix::t(Matrix::sparse.model.matrix(~ as.factor(Z) - 1))
  A1 <- Matrix::cbind2(A1_small, Matrix::Matrix(0, nrow=nrow(A1_small), ncol = aux_dim))
  n1j <- as.vector(A1_small %*% trt)
  l1 <- n1j
  u1 <- n1j
  if(verbose) message("\tx Upper and lower bounds")
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
  l2 <- sapply(Z, function(z) lowlim * n1j[z])
  u2 <- sapply(Z, function(z) uplim * n1j[z])

  if(exact_global) {
      if(verbose) message("\tx Enforce exact global balance")
      # Constrain the overall mean for individual covariates to be equal to the target
      A3 <- t(ind_covs)
      A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
      trt_sum <- colSums(ind_covs[trt == 1,, drop = F])
      l3 <- trt_sum
      u3 <- trt_sum
  } else {
      if(verbose) message("\t(SKIPPING) Enforce exact global balance")
      # skip this constraint and just make empty
      A3 <- matrix(, nrow = 0, ncol = ncol(A2))
      l3 <- numeric(0)
      u3 <- numeric(0)
  }

  if(verbose) message("\tx Fit weights to data")
  # constrain the auxiliary weights to be sqrt(P)'gamma
  sqrtP <- rbind(t(ind_covs), t(interact_covs))
  A4 <- Matrix::cbind2(sqrtP, -Matrix::Diagonal(aux_dim))
  l4 <- rep(0, aux_dim)
  u4 <- rep(0, aux_dim)

  if(verbose) message("\tx Constrain treated weights to be zero")
  # zero out treated units
  A5 <- Matrix::bdiag(lapply(trt, function(x) Matrix::Diagonal(x = x)))
  A5 <- Matrix::cbind2(A5, Matrix::Matrix(0, nrow = nrow(A5), ncol = aux_dim))
  l5 <- numeric(n)
  u5 <- numeric(n)

  if(verbose) message("\tx Combining constraints")
  A <- rbind(A1, A2, A3, A4, A5)
  l <- c(l1, l2, l3, l4, l5)
  u <- c(u1, u2, u3, u4, u5)

  return(list(A = A, l = l, u = u))
}