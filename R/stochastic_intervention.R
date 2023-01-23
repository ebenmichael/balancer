#' Balance towards a stochastic intervention
#' @param X n x d matrix of covariates
#' @param trt Vector of K-level treatment assignments 
#' @param stoch_int n x K matrix of treatment probabilities under stochastic intervention
#' @param Z Vector of group indicators with J levels
#' @inheritParams l2_balance_internal
#' @return \itemize{
#'          \item{weights }{Estimated weights as a length n vector}
#'          \item{imbalance }{Imbalance in covariates as a d X J matrix}
#'          \item{global_imbalance}{Overall imbalance in covariates, as a length d vector }}
#' @export
stochastic_int <- function(X, trt, stoch_int,
                           Z = NULL,
                           lambda = 0, lowlim = 0, uplim = 1,
                           exact_global = FALSE,
                           verbose = TRUE,
                           eps_abs = 1e-5, eps_rel = 1e-5, ...) {
  n <- nrow(X)
  K <- ncol(stoch_int)
  if(is.null(Z)) {
    Z <- rep(1, n)
  }

  # create treatment by group interaction
  trt <- as.factor(trt)
  Z <- as.factor(Z)
  trt_by_z <- interaction(trt, Z, sep = ":")
  X_trtz <- split.data.frame(X, trt_by_z)
  idxs_trtz <- split(1:nrow(X), trt_by_z)

  J <- length(X_trtz)
  aux_dim <- J * ncol(X)

  trt_vals <- sapply(strsplit(names(X_trtz), ":"), `[`, 1)
  Z_vals <- sapply(strsplit(names(X_trtz), ":"), `[`, 2)
  target_trtz <- lapply(
    1:length(X_trtz),
    function(j) {
      lev <- strsplit(levels(trt_by_z), ":")[[j]]
      trt_lev <- lev[1]
      Z_lev <- lev[2]
      idxs <- do.call(c, idxs_trtz[Z_vals == Z_lev])
      colSums(X[idxs,, drop = F] * stoch_int[idxs, trt_lev]) / n
  })
  sum_constraints <- sapply(
    1:length(X_trtz),
    function(j) {
      lev <- strsplit(levels(trt_by_z), ":")[[j]]
      trt_lev <- lev[1]
      Z_lev <- lev[2]
      idxs <- do.call(c, idxs_trtz[Z_vals == Z_lev])
      sum(stoch_int[idxs, trt_lev]) / n
  })
  target_prop_trtz <- sapply(
    1:length(X_trtz),
    function(j) {
      lev <- strsplit(levels(trt_by_z), ":")[[j]]
      trt_lev <- lev[1]
      Z_lev <- lev[2]
      idxs <- do.call(c, idxs_trtz[Z_vals == Z_lev])
      length(idxs) / n
  })
  print(do.call(c,target_trtz))
  print(sum_constraints)
  # target_prop_trtz <- sapply(idxs_trtz, length) / n
  print(target_prop_trtz)
  print(sum(sum_constraints))

  
  # Setup the components of the QP and solve
  if(verbose) message("Creating linear term vector...")
  # concenate targets for each group
  q <- - do.call(c, target_trtz)

  q <- Matrix::sparseVector(q, (n + 1):(n + aux_dim),
                            n + aux_dim)

  if(verbose) message("Creating quadratic term matrix...")
  P <- Matrix::bdiag(Matrix::Matrix(0, n, n), Matrix::Diagonal(aux_dim))

  I0 <- create_I0_matrix_multi(Xz, FALSE, n, aux_dim)
  P <- P + lambda * I0


  if(verbose) message("Creating constraint matrix...")
  constraints <- create_constraints_stochastic_int(X_trtz, target_trtz,
                                                   target_prop_trtz,
                                                   sum_constraints,
                                                   lowlim, uplim,
                                                   exact_global, verbose)


  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                              eps_rel = eps_rel,
                              eps_abs = eps_abs),
                      list(...)))

  sol <- osqp::solve_osqp(P, q, constraints$A,
                                  constraints$l, constraints$u,
                                  pars = settings)
  n_trtz <- sapply(idxs_trtz, length)
  cumsumnj <- cumsum(c(1, n_trtz))
  imbalance <- do.call(rbind, lapply(1:J,
                      function(j) {
                        wts <- sol$x[cumsumnj[j]:(cumsumnj[j + 1] - 1)]
                        c(target_trtz[[j]] - Matrix::t(X_trtz[[j]]) %*% wts)
                      }))

  if(verbose) message("Reordering weights...")
  weights <- reorder_weights(sol$x, n, NULL, trt_by_z)


  if(exact_global) {
    trt_vals <- sapply(strsplit(names(X_trtz), ":"), `[`, 1)
    global_imbal <- do.call(rbind,
        lapply(unique(trt_vals),
             function(lev) {
              Reduce(`+`,
                     lapply(1:J,
                            function(j) imbalance[j,] * target_prop_trtz[j] * (trt_vals[j] == lev)))
             }))
  } else {
    global_imbal <- NULL
  }


  # scale weights to sum to number of units in the cell
  # weights <- weights * target_prop_trtz[trt_by_z] * n / sum_constraints[trt_by_z]
  weights <- weights * n

  return(list(weights = weights,
              imbalance = imbalance,
              global_imbalance = global_imbal))
}



#' Create the constraints for QP: l <= Ax <= u
#' @param Xz length J list of n_j x d matrix of covariates for each group j
#' @param targetz length J list of d-dimensional vectors of targets for each group j
#' @param target_propz J-dimensional vector of group proportions in the target population
#' @param sum_constrains J-dimenional vector of sum constraints
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param exact_global Boolean whether to include an exact global constraint
#' @param verbose Boolean whether to display progress
#'
#' @return A, l, and u
create_constraints_stochastic_int <- function(Xz, targetz, target_propz,
                                              sum_constraints, lowlim,
                                              uplim, exact_global, verbose) {

  J <- length(Xz)

  d <- ncol(Xz[[1]])
  n <- Reduce(`+`, lapply(Xz, nrow))
  Xzt <- lapply(Xz, Matrix::t)

  # dimension of auxiliary weights
  aux_dim <- J * d



  if(verbose) message("\tx Sum to one constraint")
  # sum-to-1 constraint for each group
  A1 <- Matrix::t(Matrix::bdiag(lapply(Xz, function(x) rep(1, nrow(x)))))
  A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = aux_dim))
  l1 <- sum_constraints
  u1 <- sum_constraints
  if(verbose) message("\tx Upper and lower bounds")
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
  l2 <- rep(lowlim, n)
  u2 <- rep(uplim, n)

  if(exact_global) {
      if(verbose) message("\tx Enforce exact global balance")
      # Constrain the overall mean to be equal to the target for each treatment group
      trt_vals <- sapply(strsplit(names(Xz), ":"), `[`, 1)
      A3 <- do.call(rbind, lapply(unique(trt_vals),
             function(lev) {
              do.call(cbind,
                      lapply(1:J,
                             function(j) Xzt[[j]] * (trt_vals[j] == lev)))
             }))
      A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
      avg_target <- do.call(c,
        lapply(unique(trt_vals),
             function(lev) {
              Reduce(`+`,
                     lapply(1:J,
                            function(j) targetz[[j]] * (trt_vals[j] == lev)))
             }))
      print(avg_target)
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
