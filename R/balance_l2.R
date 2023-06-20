#' Re-weight data to a target with local and global constraints
#' @param Xz length J list of n_j x d matrix of covariates for each group j
#' @param targetz length J list of d-dimensional vectors of targets for each group j
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param exact_global Whether to enforce exact balance for overall population
#' @param target_propz J-dimensional vector of group proportions in the target population, must not be NULL if exact_global is TRUE
#' @param verbose Whether to show messages, default T
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#'
#' @return vector of weights solving balancing optimization problem
l2_balance_internal <- function(Xz, targetz,
                                lambda = 0, lowlim = 0, uplim = 1,
                                scale_sample_size = T,
                                exact_global = T, target_propz = NULL,
                                verbose = TRUE,
                                eps_abs = 1e-5, eps_rel = 1e-5, ...) {


  if(exact_global & is.null(target_propz)) {
    stop("If enforcing an exact global constraint with exact_global = T, then
         target_propz must not be NULL")
  }
  J <- length(Xz)
  aux_dim <- J * ncol(Xz[[1]])
  nz <- sapply(Xz, nrow)

  n <- sum(nz)


  # Setup the components of the QP and solve
  if(verbose) message("Creating linear term vector...")
  # concenate targets for each group
  q <- - do.call(c, targetz)

  q <- Matrix::sparseVector(q, (n + 1):(n + aux_dim),
                            n + aux_dim)

  if(verbose) message("Creating quadratic term matrix...")
  P <- Matrix::bdiag(Matrix::Matrix(0, n, n), Matrix::Diagonal(aux_dim))

  I0 <- create_I0_matrix_multi(Xz, scale_sample_size, n, aux_dim)
  P <- P + lambda * I0

  if(verbose) message("Creating constraint matrix...")
  constraints <- create_constraints_l2(Xz, targetz, target_propz, lowlim, uplim,
                                       exact_global, verbose)


  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                              eps_rel = eps_rel,
                              eps_abs = eps_abs),
                      list(...)))

  solution <- osqp::solve_osqp(P, q, constraints$A,
                                  constraints$l, constraints$u,
                                  pars = settings)

  cumsumnj <- cumsum(c(1, nz))
  imbalance <- do.call(rbind, lapply(1:J,
                      function(j) {
                        wts <- solution$x[cumsumnj[j]:(cumsumnj[j + 1] - 1)]
                        targetz[[j]] - Matrix::t(Xz[[j]]) %*% wts
                      }))
  

  if(exact_global) {
    global_imbal <- colSums(t(t(imbalance) * target_propz))
  } else {
    global_imbal <- NULL
  }
  
  # compute overall imbalance

  return(list(weights = solution$x[1:n],
              imbalance = imbalance,
              global_imbalance = global_imbal
              ))

}


#' Create the constraints for QP: l <= Ax <= u
#' @param Xz length J list of n_j x d matrix of covariates for each group j
#' @param targetz length J list of d-dimensional vectors of targets for each group j
#' @param target_propz J-dimensional vector of group proportions in the target population
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param exact_global Boolean whether to include an exact global constraint
#' @param verbose Boolean whether to display progress
#'
#' @return A, l, and u
create_constraints_l2 <- function(Xz, targetz, target_propz, lowlim, uplim,
                                 exact_global, verbose) {

  J <- length(Xz)

  d <- ncol(Xz[[1]])
  n <- Reduce(`+`, lapply(Xz, nrow))
  Xzt <- lapply(Xz, Matrix::t)

  # dimension of auxiliary weights
  aux_dim <- J * d



  if(verbose) message("\tx Sum to one constraint")
  # sum-to-target proportions constraint for each group
  A1 <- Matrix::t(Matrix::bdiag(lapply(Xz, function(x) rep(1, nrow(x)))))
  A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = aux_dim))
  l1 <- rep(1, J)
  u1 <- rep(1, J)
  if(verbose) message("\tx Upper and lower bounds")
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
  l2 <- rep(lowlim, n)
  u2 <- rep(uplim, n)

  if(exact_global) {
      if(verbose) message("\tx Enforce exact global balance")
      # Constrain the overall mean to be equal to the target
      A3 <- do.call(cbind, lapply(1:J, function(j) Xzt[[j]] * target_propz[j]))
      A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
      avg_target <- Reduce(`+`,
          lapply(1:J, function(j) target_propz[j] * targetz[[j]])
        )
      l3 <- avg_target
      u3 <- avg_target
      
  } else {
      # if(verbose) message("\t(SKIPPING) Enforce exact global balance")
      # skip this constraint and just make empty
      A3 <- matrix(, nrow = 0, ncol = ncol(A2))
      l3 <- numeric(0)
      u3 <- numeric(0)
  }

    if(verbose) message("\tx Fit weights to data")
    # constrain the auxiliary weights to be sqrt(P)'gamma
    sqrtP <- Matrix::bdiag(Xzt)
    A4 <- Matrix::cbind2(sqrtP, -Matrix::Diagonal(aux_dim))
    l4 <- rep(0, aux_dim)
    u4 <- rep(0, aux_dim)

    if(verbose) message("\tx Combining constraints")
    A <- rbind(A1, A2, A3, A4)
    l <- c(l1, l2, l3, l4)
    u <- c(u1, u2, u3, u4)

    return(list(A = A, l = l, u = u))
}


### Next are a series of special cases of balance_l2_internal


reorder_weights <- function(sol, n, trtz, Z) {

  if(is.null(trtz)) {
    trtz <- split(numeric(n), Z)
  }

  # convert weights into a matrix
  J <- length(trtz)
  nz0 <- sapply(1:J, function(j) sum(1 - trtz[[j]]))
  nz <- sapply(trtz, length)
  weights <- numeric(n)
  idxs <- split(1:length(Z), Z)

  cumsumnj <- cumsum(c(1, nz0))
  for(j in 1:J) {
    weightsj <- numeric(nz[j])
    weightsj[trtz[[j]] == 0] <- sol[cumsumnj[j]:(cumsumnj[j + 1] - 1)]
    weights[idxs[[j]]] <- weightsj
  }

  return(weights)
}







# #' Re-weight control sub-groups to treated sub-group means
# #' @param X n x d matrix of covariates
# #' @param trt Vector of treatment assignments
# #' @param Z Vector of group indicators with J levels
# #' @inheritParams l2_balance_internal
# #'
# #' @return \itemize{
# #'          \item{weights }{Estimated weights as a length n vector}
# #'          \item{imbalance }{Imbalance in covariates as a d X J matrix}
# #'          \item{global_imbalance}{Overall imbalance in covariates, as a length d vector }}
# #' @export
# multilevel_qp2 <- function(X, trt, Z, lambda = 0, lowlim = 0, uplim = 1,
#                         scale_sample_size = T, exact_global = T,
#                         verbose = TRUE,
#                         eps_abs = 1e-5, eps_rel = 1e-5, ...) {


  
#   # convert X to a matrix
#   X <- as.matrix(X)

#   # split data and treatment by factor
#   Z_factor <- as.factor(Z)
#   Xz <- split.data.frame(X, Z_factor)
#   trtz <- split(trt, Z)
#   J <- length(Xz)
  

#   check_data_multi(X, trt, Z, Xz, lambda, lowlim, uplim)

#   # create targets  and target proportions as treatment averages
#   targetz <- lapply(1:J, function(j) colMeans(Xz[[j]][trtz[[j]] == 1,, drop = F]))
#   n1z <- sapply(trtz, sum)
#   target_propz <- n1z / sum(n1z)

#   # get control units in each group
#   Xz_ctrl <- lapply(1:J, function(j) Xz[[j]][trtz[[j]] == 0,, drop = F])

  
#   sol <- l2_balance_internal(Xz_ctrl, targetz, lambda, lowlim, uplim,
#                              scale_sample_size,
#                              exact_global, target_propz,
#                              verbose, eps_abs, eps_rel, ...)

#   if(verbose) message("Reordering weights...")
#   weights <- reorder_weights(sol$weights, n, trtz, Z_factor)
#   # scale weights to sum to n1j
#   weights <- weights * n1z[Z]

#   return(list(weights = weights,
#               imbalance = sol$imbalance,
#               global_imbalance = sol$global_imbal))
# }
