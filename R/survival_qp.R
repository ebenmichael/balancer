################################################################################
## Balancing weights for survival analysis
################################################################################

#' Re-weight control sub-groups to treated sub-group means
#' @param B_X n x k basis matrix of covariates 
#' @param trt Vector of treatment assignments
#' @param times Vector of event/censoring times (see events for if event or censoring time)
#' @param events Vector of boolean censoring indicators (whether individual was censored)
#' @param t Time
#' @param lambda Regularization hyperparameter, default 0
#' @param lowlim Lower limit on weights, default 1
#' @param uplim Upper limit on weights, default NULL
#' @param verbose Whether to show messages, default T
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#'
#' @return \itemize{
#'          \item{weights}{Estimated weights as a length n vector}
#'          \item{imbalance1}{Imbalance in covariates for treated in vector}
#'          \item(imbalance0){Imbalance in covariates for control in vector}}
#' @export
survival_qp <- function(B_X, trt, times, events, t, lambda = 0, lowlim = 1, uplim = NULL,
                        verbose = TRUE, eps_abs = 1e-5, eps_rel = 1e-5, ...) {
  # Convert X to a matrix
  B_X <- as.matrix(B_X)
  
  # Simple data checks
  check_data_surv(B_X, trt, times, events, t, lambda, lowlim, uplim)
  
  n <- nrow(B_X)
  # If no uplim set (default NULL), set uplim to be number of individuals n
  if (is.null(uplim)) {
    uplim <- n
  }
  
  # Calculate mean basis vector from basis mx on covariates 
  Bbar_X <- colMeans(B_X)
  
  # Setup the components of the QP and solve
  if(verbose) message("Creating linear term vector...")
  q <- create_q_vector_surv(n, trt, B_X, Bbar_X)
  
  if(verbose) message("Creating quadratic term matrix...")
  P <- create_P_matrix_surv(n, trt, B_X)
  
  I0 <- create_I0_matrix_surv(n)
  
  P <- P + lambda * I0
  
  if(verbose) message("Creating constraint matrix...")
  constraints <- create_constraints(n, trt, times, events, t, lowlim, uplim, verbose)
  
  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                             eps_rel = eps_rel,
                             eps_abs = eps_abs),
                        list(...)))
  
  solution <- osqp::solve_osqp(P, q, constraints$A,
                               constraints$l, constraints$u,
                               settings)
  
  weights <- solution$x

  # Return vectors of imbalances for treated and control 
  noncens_t <- ((times >= t) & (events == TRUE)) | (events == FALSE)
  B0_X <- B_X*(trt == 0)
  B1_X <- B_X*(trt == 1)
  
  imbalance <- rowSums(sweep((weights * noncens_t)*B0_X, 2, Bbar_X)^2) + rowSums(sweep((weights * noncens_t)*B1_X, 2, Bbar_X)^2)
  
  # imbalance <- rowSums(sweep((weights * noncens_t)* B_X, 2, Bbar_X)^2)
  
  return(list(weights = weights,
              imbalance1 = imbalance[trt == 1],
              imbalance0 = imbalance[trt == 0]))
}

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param B_X n x k basis matrix of covariates 
#' @param Bbar_X Vector of population means to re-weight to
#'
#' @return q vector
create_q_vector_surv <- function(n, trt, B_X, Bbar_X) {
  B0_X <- B_X*(trt == 0)
  B1_X <- B_X*(trt == 1)

  q <- -(2/n) * Bbar_X %*% t(B0_X + B1_X) 
  return(q)
}

#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#' @param B_X n x k basis matrix of covariates 
#'
#' @return P matrix
create_P_matrix_surv <- function(n, trt, B_X) {
  B0_X <- B_X*(trt == 0)
  B0B0 <- B0_X %*% t(B0_X)
  
  B1_X <- B_X*(trt == 1)
  B1B1 <- B1_X %*% t(B1_X)
  
  P <- (1/n^2) * (B0B0 + B1B1)
  return(P)
}

#' Create diagonal regularization matrix
#' @param n Total number of units
#' @param lambda Regularization hyperparameter, default 0
#' 
#' @return I0 matrix
create_I0_matrix_surv <- function(n, lambda) {
  
  I0 <- diag(n)
  return(I0)
}

#' Create the constraints for QP: l <= Ax <= u
#' @param n Number of individuals
#' @param trt Vector of treatment assignments
#' @param times Vector of event/censoring times (see events for if event or censoring time)
#' @param events Vector of boolean censoring indicators (whether individual was censored)
#' @param t Time
#' @param lowlim Lower limit on weights, default 1
#' @param uplim Upper limit on weights, default NULL
#' @param verbose Whether to show messages, default T
#' 
#' @return A, l, and u
create_constraints <- function(n, trt, times, events, t, lowlim = 1, uplim = NULL, verbose = TRUE) {
  
  # Vector of booleans indicating if indiv is not censored at time t
  noncens_t <- ((times >= t) & (events == TRUE)) | (events == FALSE)
  
  if(verbose) message("\tx Summed weight constraint for non-censored treatment group")
  A1 <- as.integer(noncens_t & (trt == 1))
  l1 <- n
  u1 <- n
  
  if(verbose) message("\tx Summed weight constraint for non-censored control group")
  A2 <- as.integer(noncens_t & (trt == 0))
  l2 <- n
  u2 <- n
  
  if(verbose) message("\tx Restrict weights for non-censored individuals (between 1 and n)")
  A3 <- diag(as.integer(noncens_t))
  l3 <- rep(1, n)
  u3 <- rep(n, n)
  
  if(verbose) message("\tx Restrict weights for censored individuals (equal to 0)")
  A4 <- diag(as.integer(!noncens_t))
  l4 <- rep(0, length(noncens_t))
  u4 <- rep(0, length(noncens_t))
  
  if(verbose) message("\tx Combining constraints")
  A <- rbind(A1, A2, A3, A4)
  l <- c(l1, l2, l3, l4)
  u <- c(u1, u2, u3, u4)
  
  if(verbose) message("\tx Removing null constraints")
  nonNullConstrIdxs <- which(rowSums(A) != 0)
  A <- A[nonNullConstrIdxs, ]
  l <- l[nonNullConstrIdxs]
  u <- u[nonNullConstrIdxs]
  
  return(list(A = A, l = l, u = u))
}

#' Check that data is in right shape, hyparameters are feasible, and basis fns are valid
#' @param B_X n x k basis matrix of covariates (X)
#' @param trt Vector of treatment assignments
#' @param times Vector of event/censoring times (see events for if event or censoring time)
#' @param events Vector of boolean censoring indicators (whether individual was censored)
#' @param t Time
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 1
#' @param uplim Upper limit on weights, default NULL
check_data_surv <- function(B_X, trt, times, events, t, lambda, lowlim = 1, uplim = NULL) {
  # NA checks
  if(any(is.na(B_X))) {
    stop("Basis matrix B_X contains NA values")
  }
  
  if(any(! trt %in% c(0,1))) {
    stop("Treatment must be (0,1)")
  }
  
  if(any(! events %in% c(TRUE,FALSE))) {
    stop("Censoring indicators must be (TRUE,FALSE)")
  }
  
  if(any(is.na(times))) {
    stop("NA values in times")
  }
  
  # Group size check
  # Vector of booleans indicating if indiv is not censored at time t
  # TRUE: indiv is not censored at time t, FALSE: indiv is censored at time t
  noncens_t <- ((times >= t) & (events == TRUE)) | (events == FALSE)
  if(sum(noncens_t & (trt == 0)) == 0) {
    stop("No non-censored individuals in treatment group at time ", t)
  }
  if(sum(noncens_t & (trt == 1)) == 0) {
    stop("No non-censored individuals in control group at time ", t)
  }
  
  #dimension checks
  n <- nrow(B_X)
  # k <- ncol(B_X)
  
  if (length(trt) != n) {
    stop("Number of rows in basis matrix B_X (",
         n, ") does not equal the length of trt (",
         length(trt), ")")
  }
  
  # hyerparameters are feasible
  if(lambda < 0) {
    stop("lambda must be >= 0")
  }
  if(lowlim > uplim) {
    stop("Lower threshold must be lower than upper threshold")
  }
  # if(lowlim < 0) {
  #   stop("Lower threshold must be at least 0")
  # }
  # if(uplim > n) {
  #   stop("Upper threshold must be at most n (",
  #        n, ")")
  # }
}