################################################################################
## Wrapper to standardize to target means
################################################################################

#' Re-weight groups to target population means
#' @param X0 n x d0 matrix of (transformed) covariates defining the mean control response function
#' @param Xtau n x dtau matrix of (transformed) covariates defining the mean treatment effect function
#' @param target Vector of population means to re-weight to
#' @param S Numeric vector of site indicators with J levels
#' @param Z Numeric vector of treatment indicators with 2 levels
#' @param pscores Numeric vector of propensity scores
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default nrow(X0)
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default F
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_program Whether to return the objective matrix and vector and constraints, default T
#' @param init_uniform Wheter to initialize solver with uniform weights, default F
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#'
#' @return \itemize{
#'          \item{weights }{Estimated primal weights as an n x J matrix}
#'          \item{imbalance_0 }{Imbalance in covariates defining mean control response function as a d0 x J matrix}
#'          \item{imbalance_tau }{Imbalance in covariates defining mean treatment effect function as a dtau x J matrix}
#'          \item{data_out }{List containing elements of QP min 0.5 x'Px + q'x st l <= Ax <= u \itemize{
#'                  \item{P, q}{}
#'                  \item{constraints }{A, l , u}
#'}}}
#' @export
standardize_treatment <- function(X0, Xtau, target, S, Z, pscores, lambda = 0,
                                  lowlim = 0, uplim = nrow(X0), scale_sample_size = F,
                                  data_in = NULL, verbose = TRUE,
                                  return_program = TRUE,
                                  init_uniform = F, eps_abs = 1e-5,
                                  eps_rel = 1e-5, ...) {

  # ensure that covariate matrices are matrices and get total number of units
  X0 <- as.matrix(X0)
  Xtau <- as.matrix(Xtau)
  n <- nrow(X0)

  # convert site indicators to factor, get index of sites, and compute # sites
  S_factor <- as.factor(S)
  unique_S <- levels(S_factor)
  J <- length(unique_S)

  # split covariate matrix by site and compute number of units in each site
  X0s <- split.data.frame(X0, S_factor)
  Xtaus <- split.data.frame(Xtau, S_factor)
  nj <- sapply(X0s, nrow)

  # get row IDs of units by site
  idxs <- split(1:nrow(X0), S_factor)

  # ensure that target is a vector
  target <- c(target)

  # check arguments for issues
  check_args_treatment(X0, Xtau, target, S, Z, pscores, X0s, Xtaus,
                       as.numeric(nj), lambda, lowlim, uplim, data_in)

  # create propensity score multipliers and split by site
  pro_trt <- Z / (pscores * nj[S_factor])
  pro_trt_split <- split(pro_trt, S_factor)
  pro_ctr <- (1 - Z) / ((1 - pscores) * nj[S_factor])
  pro_ctr_split <- split(pro_ctr, S_factor)

  # define number of auxiliary weights
  aux_dim <- (J * ncol(X0)) + (J * ncol(Xtau))

  # construct linear term vector
  if(verbose) message("Creating linear term vector...")
  if(is.null(data_in$q)) {
    q <- create_q_vector_treatment(Xtaus, pro_trt_split, target, aux_dim)
  } else {
    q <- data_in$q
  }

  # construct quadratic term matrix
  if(verbose) message("Creating quadratic term matrix...")
  if(is.null(data_in$P)) {
    P <- create_P_matrix_treatment(n, aux_dim)
  } else {
    P <- data_in$P
  }
  I0 <- create_I0_matrix_treatment(pro_trt_split, pro_ctr_split,
                                   scale_sample_size, nj, n, aux_dim)
  P <- P + lambda * I0

  # construct constraint matrix
  if(verbose) message("Creating constraint matrix...")
  if(is.null(data_in$constraints)) {
    constraints <- create_constraints_treatment(X0s, Xtaus, target, Z, S_factor,
                                                pro_trt_split, pro_ctr_split,
                                                lowlim, uplim, verbose)
  } else {
    constraints <- data_in$constraints
    constraints$l[(J + 1):(J + n)] <- lowlim
    constraints$u[(J + 1):(J + n)] <- uplim
  }

  # set optimization settings
  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                             eps_rel = eps_rel,
                             eps_abs = eps_abs),
                        list(...)))

  # solve optimization problem (possibly with uniform weights)
  if(init_uniform) {
    if(verbose) message("Initializing with uniform weights")
    unifw <- get_uniform_weights_treatment(X0s, Xtaus)
    obj <- osqp::osqp(P, q, constraints$A,
                      constraints$l, constraints$u, pars = settings)
    obj$WarmStart(x = unifw)
    solution <- obj$Solve()
  } else {
    solution <- osqp::solve_osqp(P, q, constraints$A,
                                 constraints$l, constraints$u,
                                 pars = settings)
  }

  # convert weights into a matrix
  if(verbose) message("Reordering weights...")
  weights <- matrix(0, ncol = J, nrow = n)
  cumsumnj <- cumsum(c(1, nj))
  for(j in 1:J) {
    weights[idxs[[j]], j] <- solution$x[cumsumnj[j]:(cumsumnj[j + 1] - 1)]
  }

  # compute imbalance matrix
  imbalancetau <- as.matrix(target - t(Xtau) %*% (weights * pro_trt))
  imbalance0 <- as.matrix(t(X0) %*% (weights * pro_trt) -
                          t(X0) %*% (weights * pro_ctr))

  # collapse weight matrix to vector
  weights <- rowSums(weights)

  # package program components if requested by user
  if(return_program) {
    program <- list(P = P  - lambda * I0,
                     q = q, constraints = constraints)
  } else {
    program <- NULL
  }

  # return output
  return(list(weights = weights, imbalance_0 = imbalance0,
              imbalance_tau = imbalancetau, program = program))

}


#' Create diagonal regularization matrix
#' @param pro_trt_split List of J vectors of treatment propensity multipliers
#' @param pro_ctr_split List of J vectors of control propensity multipliers
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param nj Number of units in each group
#' @param n Total number of units
#' @param aux_dim Dimension of auxiliary weights
create_I0_matrix_treatment <- function(pro_trt_split, pro_ctr_split, scale_sample_size, nj, n, aux_dim) {

  if(scale_sample_size) {
    # diagonal matrix of inverse propensity scores scaled by group sample size
    I0 <- Matrix::Diagonal(n, rep(nj, nj) * (unlist(pro_trt_split) + unlist(pro_ctr_split)))
  } else {
    # diagonal matrix of inverse propensity scores not scaled by group sample size
    I0 <- Matrix::Diagonal(n, unlist(pro_trt_split) + unlist(pro_ctr_split))
  }
  # add placeholder 0s for auxiliary weight entries
  I0 <- Matrix::bdiag(I0, Matrix::Diagonal(aux_dim, 0))
  return(I0)
}

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param Xtaus List of J n x dtau matrices of covariates split by group
#' @param pro_trt_split List of J vectors of treatment propensity multipliers
#' @param target Vector of population means to re-weight to
#' @param aux_dim Dimension of auxiliary weights
#'
#' @return q vector
create_q_vector_treatment <- function(Xtaus, pro_trt_split, target, aux_dim) {
  q <- -c(do.call(rbind, Xtaus) %*% target) * unlist(pro_trt_split)
  q <- Matrix::sparseVector(q, 1:length(q),
                            length(q) + aux_dim)
  return(q)
}


#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#' @param n Total number of units
#' @param aux_dim Dimension of auxiliary weights
#'
#' @return P matrix
create_P_matrix_treatment <- function(n, aux_dim) {
  return(Matrix::bdiag(Matrix::Diagonal(n, 0), Matrix::Diagonal(aux_dim, 1)))
}

#' Get a set of uniform weights for initialization
#' @param X0s List of J n x d0 matrices of covariates split by group
#' @param Xtaus List of J n x dtau matrices of covariates split by group
#'
get_uniform_weights_treatment <- function(X0s, Xtaus) {

  # uniform weights for each group
  uniw <- do.call(c, lapply(X0s, function(x) rep(1 / nrow(x), nrow(x))))

  # transformed auxiliary uniform weights
  sqrtP0 <- Matrix::bdiag(lapply(X0s, t))
  aux0_uniw <- as.numeric(sqrtP0 %*% uniw)

  sqrtPtau <- Matrix::bdiag(lapply(Xtaus, t))
  auxtau_uniw <- as.numeric(sqrtPtau %*% uniw)
  return(c(uniw, c(aux0_uniw, auxtau_uniw)))
}

#' Create the constraints for QP: l <= Ax <= u
#' @param X0s List of J n x d0 matrices of covariates split by group
#' @param Xtaus List of J n x dtau matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators
#' @param S_factor Vector of site indicators
#' @param pro_trt_split List of J vectors of treatment propensity multipliers
#' @param pro_ctr_split List of J vectors of control propensity multipliers
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param verbose Boolean indicating whether to print progress
#'
#' @return A, l, and u
create_constraints_treatment <- function(X0s, Xtaus, target, Z, S_factor,
                               pro_trt_split, pro_ctr_split,
                               lowlim, uplim, verbose) {

  J <- length(X0s)
  nj <- as.numeric(sapply(X0s, nrow))
  d0 <- ncol(X0s[[1]])
  dtau <- ncol(Xtaus[[1]])
  n <- sum(nj)
  X0st <- lapply(X0s, t)
  Xtaust <- lapply(Xtaus, t)

  # compute number of treated and control units within each site
  n1j <- sapply(unique(S_factor), FUN = function(j) sum(Z[S_factor == j]))
  n0j <- sapply(unique(S_factor), FUN = function(j) sum(1 - (Z[S_factor == j])))

  # dimension of auxiliary weights
  aux_dim <- (J * d0) + (J * dtau)

  if(verbose) message("\tx Sum to one constraint")
  # sum-to-one constraint for each group
  trt_mat <- Matrix::t(Matrix::bdiag(split(Z, S_factor)))
  ctr_mat <- Matrix::t(Matrix::bdiag(split(1 - Z, S_factor)))
  A1 <- Matrix::rbind2(trt_mat, ctr_mat)

  rm(trt_mat)
  rm(ctr_mat)

  A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = aux_dim))
  l1 <- c(n1j, n0j)
  u1 <- c(n1j, n0j)

  if(verbose) message("\tx Upper and lower bounds")
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
  l2 <- rep(lowlim, n)
  u2 <- rep(uplim, n)

  if(verbose) message("\tx Fit weights to data")
  # constrain the auxiliary weights to be sqrt(P)'gamma
  sqrtP1 <- Matrix::t(Matrix::t(Matrix::bdiag(X0st)) * (unlist(pro_trt_split) - unlist(pro_ctr_split)))
  sqrtP2 <- Matrix::t(Matrix::t(Matrix::bdiag(Xtaust)) * unlist(pro_trt_split))
  A3 <- Matrix::cbind2(Matrix::rbind2(sqrtP1, sqrtP2), -Matrix::Diagonal(aux_dim))
  l3 <- rep(0, aux_dim)
  u3 <- rep(0, aux_dim)

  if(verbose) message("\tx Combining constraints")
  A <- rbind(A1, A2, A3)
  l <- c(l1, l2, l3)
  u <- c(u1, u2, u3)

  return(list(A = A, l = l, u = u))
}



#' Check that data is in right shape and hyparameters are feasible
#' @param X0 n x d0 matrix of covariates
#' @param Xtau n x dtau matrix of covariates
#' @param target Vector of population means to re-weight to
#' @param S Vector of group indicators with J levels
#' @param Z Vector of treatment indicators
#' @param pscores Vector of propensity scores
#' @param X0s List of J n x d0 matrices of covariates split by group
#' @param Xtaus List of J n x dtau matrices of covariates split by group
#' @param nj Number of units in jth site
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
check_args_treatment <- function(X0, Xtau, target, S, Z, pscores, X0s, Xtaus, nj, lambda, lowlim, uplim, data_in) {

  # NA checks
  if(any(is.na(X0))) {
    stop("Covariate matrix X0 contains NA values.")
  }

  if(any(is.na(S))) {
    stop("Grouping vector S contains NA values.")
  }

  if(any(is.na(Z))) {
    stop("Treatment vector Z contains NA values.")
  }

  if(any(is.na(pscores))) {
    stop("Propensity score vector pscores contains NA values.")
  }

  if(any(is.na(target))) {
    stop("Target vector contains NA values.")
  }

  # dimension checks
  n <- nrow(X0)
  d0 <- ncol(X0)
  dtau <- ncol(Xtau)
  J <- length(X0s)
  aux_dim <- (J * d0) + (J * dtau)

  if(length(S) != n) {
    stop("The number of rows in covariate matrix X (", n,
         ") does not equal the length of group vector S (",
         length(S), ").")
  }

  if(length(Z) != n) {
    stop("The number of rows in covariate matrix X (", n,
         ") does not equal the length of treatment vector Z (",
         length(Z), ").")
  }

  if(length(pscores) != n) {
    stop("The number of rows in covariate matrix X (", n,
         ") does not equal the length of propensity score vector pscores (",
         length(pscores), ").")
  }

  if(sum(nj) != n) {
    stop("Implied number of weights (", sum(nj),
         ") does not equal number of units (", n, ").")
  }

  if(length(target) != dtau) {
    stop("Target dimension (", length(target),
         ") is not equal to data dimension (", dtau, ").")
  }

  if(!is.null(data_in$q)) {
    if(length(data_in$q) != n + aux_dim) {
      stop("data_in$q vectors should have dimension ", n + aux_dim)
    }
  }

  if(!is.null(data_in$P)) {
    if(dim(data_in$P)[1] != dim(data_in$P)[2]) {
      stop("data_in$P matrix must be square")
    }
    if(dim(data_in$P)[1] != n + aux_dim) {
      stop("data_in$P should have ", n + aux_dim,
           " rows and columns")
    }
  }

  if(!is.null(data_in$constraints)) {
    if(length(data_in$constraints$l) != length(data_in$constraints$u)) {
      stop("data_in$constraints$l and data_in$constraints$u",
           " must have the same dimension")
    }
    if(length(data_in$constraints$l) != 2 * J + n + d0 + dtau + aux_dim) {
      stop("data_in$constraints$l must have dimension ",
           2 * J + d0 + dtau + n + aux_dim)
    }
    if(nrow(data_in$constraints$A) != length(data_in$constraints$l)) {
      stop("The number of rows in data_in$constraints$A must be ",
           "the same as the dimension of data_in$constraints$l")
    }

    if(ncol(data_in$constraints$A) != n + aux_dim) {
      stop("The number of columns in data_in$constraints$A must be ",
           n + aux_dim)
    }
  }

  # hyerparameters are feasible
  if(lambda < 0) {
    stop("lambda must be >= 0")
  }
  if(lowlim > uplim) {
    stop("Lower threshold must be lower than upper threshold")
  }
  if(lowlim > 1/max(nj)) {
    stop("Lower threshold must be lower than 1 / size of largest group")
  }
  if(uplim < 1 / min(nj)) {
    stop("Upper threshold must be higher than 1 / size of smallest group")
  }
}
