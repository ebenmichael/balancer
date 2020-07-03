################################################################################
## Wrapper to standardize to target means
################################################################################

#' Re-weight groups to target population means
#' @param X0 n x d0 matrix of covariates defining the mean control response function
#' @param Xtau n x dtau matrix of covariates defining the mean treatment effect function
#' @param target Vector of population means to re-weight to
#' @param S Numeric vector of site indicators with J levels
#' @param Z Numeric vector of treatment indicators with 2 levels
#' @param pscores Numeric vector of propensity scores
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_program Whether to return the objective matrix and vector and constraints, default T
#' @param exact_global Whether to enforce exact balance for overall population
#' @param init_uniform Wheter to initialize solver with uniform weights
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#'
#' @return \itemize{
#'          \item{weights }{Estimated primal weights as an n x J matrix}
#'          \item{imbalance }{Imbalance in covariates as a d X J matrix}
#'          \item{data_out }{List containing elements of QP min 0.5 x'Px + q'x st l <= Ax <= u \itemize{
#'                  \item{P, q}{}
#'                  \item{constraints }{A, l , u}
#'}}}
#' @export
standardize_treatment <- function(X0, Xtau, target, S, Z, pscores, lambda = 0,
                                  lowlim = 0, uplim = 1, scale_sample_size = T,
                                  data_in = NULL, verbose = TRUE,
                                  return_program = TRUE, exact_global = F,
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
  check_args(X0, Xtau, target, S, Z, pscores, X0s, Xtaus, as.numeric(nj),
             lambda, lowlim, uplim, data_in)
  
  # create propensity score multipliers and split by site
  pro_trt_split <- split(Z / (pscores * nj[S_factor]), S_factor)
  pro_ctr_split <- split((1 - Z) / ((1 - pscores) * nj[S_factor]), S_factor)
  
  # define number of auxiliary weights
  aux_dim <- (J * ncol(X0)) + (J * ncol(Xtau))
  
  # construct linear term vector
  if(verbose) message("Creating linear term vector...")
  if(is.null(data_in$q)) {
    q <- create_q_vector(Xtaus, pro_trt_split, target, aux_dim)
  } else {
    q <- data_in$q
  }
  
  # construct quadratic term matrix
  if(verbose) message("Creating quadratic term matrix...")
  if(is.null(data_in$P)) {
    P <- create_P_matrix(n, aux_dim)
  } else {
    P <- data_in$P
  }
  I0 <- create_I0_matrix(pro_trt_split, pro_ctr_split,
                         scale_sample_size, nj, n, aux_dim)
  P <- P + lambda * I0
  
  # construct constraint matrix
  if(verbose) message("Creating constraint matrix...")
  if(is.null(data_in$constraints)) {
    constraints <- create_constraints(X0s, Xtaus, target, Z, S_factor,
                                      pro_trt_split, pro_ctr_split,
                                      lowlim, uplim, exact_global, verbose)
  } else {
    constraints <- data_in$constraints
    constraints$l[(J + 1):(J + n)] <- lowlim
    constraints$u[(J + 1):(J + n)] <- uplim
  }
  
  # set optimization settings
  if(verbose) {
    settings <- osqp::osqpSettings(verbose = TRUE, eps_abs = 1e-6,
                                   eps_rel = 1e-6, max_iter = 5000)
  } else {
    settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-6,
                                   eps_rel = 1e-6, max_iter = 5000)
  }
  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                             eps_rel = eps_rel,
                             eps_abs = eps_abs),
                        list(...)))

  # solve optimization problem (possibly with uniform weights)
  if(init_uniform) {
    if(verbose) message("Initializing with uniform weights")
    unifw <- get_uniform_weights(Xs)
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
  imbalancetau <- as.matrix(target - t(Xtau) %*% (weights * unlist(pro_trt_split)))
  imbalance0 <- as.matrix(t(X0) %*% (weights * unlist(pro_trt_split)) -
                          t(X0) %*% (weights * unlist(pro_ctr_split)))
  
  # package program components if requested by user
  if(return_program) {
    program <- list(P = P  - lambda * I0,
                     q = q, constraints = constraints)
  } else {
    program <- NULL
  }
  
  # return output
  return(list(weights = weights, imbalance_tau = imbalancetau,
              imbalance_0 = imbalance0, program = program))
  
}


#' Create diagonal regularization matrix
#' @param pro_trt_split List of J vectors of treatment propensity multipliers
#' @param pro_ctr_split List of J vectors of control propensity multipliers
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param nj Number of units in each group
#' @param n Total number of units
#' @param aux_dim Dimension of auxiliary weights
create_I0_matrix <- function(pro_trt_split, pro_ctr_split, scale_sample_size, nj, n, aux_dim) {
  
  if(scale_sample_size) {
    # diagonal matrix of inverse propensity scores scaled by group sample size
    I0 <- Matrix::Diagonal(n, rep(nj ^ 2, nj) * (unlist(pro_trt_split) + unlist(pro_ctr_split)))
  } else {
    # diagonal matrix of inverse propensity scores not scaled by group sample size
    I0 <- Matrix::Diagonal(n, rep(nj, nj) * (unlist(pro_trt_split) + unlist(pro_ctr_split)))
  }
  # add placeholder 0s for auxiliary weight entries
  I0 <- Matrix::bdiag(I0, Matrix::Diagonal(aux_dim, 0))
  return(I0)
}

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param Xs List of J n x d matrices of covariates split by group
#' @param pro_trt_split List of J vectors of treatment propensity multipliers
#' @param target Vector of population means to re-weight to
#' @param aux_dim Dimension of auxiliary weights
#'
#' @return q vector
create_q_vector <- function(Xtaus, pro_trt_split, target, aux_dim) {
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
create_P_matrix <- function(n, aux_dim) {
  return(Matrix::bdiag(Matrix::Diagonal(n, 0), Matrix::Diagonal(aux_dim, 1)))
}

#' Get a set of uniform weights for initialization
#' @param Xs List of J n x d matrices of covariates split by group
#' 
get_uniform_weights <- function(Xs) {
  
  # uniform weights for each group
  uniw <- do.call(c, lapply(Xs, function(x) rep(1 / nrow(x), nrow(x))))
  
  # transformed auxiliary uniform weights
  sqrtP <- Matrix::bdiag(lapply(Xs, t))
  aux_uniw <- as.numeric(sqrtP %*% uniw)
  return(c(uniw, aux_uniw))
}

#' Create the constraints for QP: l <= Ax <= u
#' @param X0s List of J n x d0 matrices of covariates split by group
#' @param Xtaus List of J n x dtau matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators
#' @param S_factor Vector of site indicators
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param exact_global Boolean indicating whether to enforce global constraint
#' @param verbose Boolean indicating whether to print progress
#'
#' @return A, l, and u
create_constraints <- function(X0s, Xtaus, target, Z, S_factor,
                               pro_trt_split, pro_ctr_split,
                               lowlim, uplim, exact_global, verbose) {
  
  J <- length(X0s)
  nj <- as.numeric(sapply(X0s, nrow))
  d0 <- ncol(X0s[[1]])
  dtau <- ncol(Xtaus[[1]])
  cumsum_nj <- cumsum(c(1, nj))
  n <- sum(nj)
  X0st <- lapply(X0s, t)
  Xtaust <- lapply(Xtaus, t)
  
  # compute number of treated and control units within each site
  n1j <- sapply(1:J, FUN = function(j) sum(Z[S == j]))
  n0j <- sapply(1:J, FUN = function(j) sum(1 - (Z[S == j])))
  
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
  
  if(exact_global) {
    if(verbose) message("\tx Mantain overall population mean")
    # require that the weighted average of the treated units (over all sites)
    # equals the weighted average of the control units (over all sites)
    response_mat <- t(
                      t(
                        do.call(cbind, X0st)) * ((unlist(pro_trt_split) * rep(n1j / sum(n1j), nj)) -
                                                 (unlist(pro_ctr_split) * rep(n0j / sum(n0j), nj))))
    treat_mat <- t(t(do.call(cbind, Xtaust)) * rep(n1j, nj) * unlist(pro_trt_split))
    
    A3 <- Matrix::rbind2(response_mat, treat_mat)
    
    rm(response_mat)
    rm(treat_mat)
    
    A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
    
    l3 <- c(rep(0, d0), sum(n1j) * target)
    u3 <- c(rep(0, d0), sum(n1j) * target)
  } else {
    if(verbose) message("\t(SKIPPING) Mantain overall population mean")
    # skip this constraint and just make empty
    A3 <- matrix(, nrow = 0, ncol = ncol(A2))
    l3 <- numeric(0)
    u3 <- numeric(0)
  }
  
  if(verbose) message("\tx Fit weights to data")
  # constrain the auxiliary weights to be sqrt(P)'gamma
  sqrtP1 <- Matrix::t(Matrix::t(Matrix::bdiag(X0st)) * (unlist(pro_trt_split) - unlist(pro_ctr_split)))
  sqrtP2 <- Matrix::t(Matrix::t(Matrix::bdiag(Xtaust)) * unlist(pro_trt_split))
  A4 <- Matrix::cbind2(Matrix::rbind2(sqrtP1, sqrtP2), -Matrix::Diagonal(aux_dim))
  l4 <- rep(0, aux_dim)
  u4 <- rep(0, aux_dim)
  
  if(verbose) message("\tx Combining constraints")
  A <- rbind(A1, A2, A3, A4)
  l <- c(l1, l2, l3, l4)
  u <- c(u1, u2, u3, u4)
  
  return(list(A = A, l = l, u = u))
}



#' Check that data is in right shape and hyparameters are feasible
#' @param X n x d matrix of covariates
#' @param target Vector of population means to re-weight to
#' @param S Vector of group indicators with J levels
#' @param Z Vector of treatment indicators
#' @param pscores Vector of propensity scores
#' @param Xs List of J n x d matrices of covariates split by group
#' @param nj Number of units in jth site
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_data Whether to return the objective matrix and vector and constraints, default T
check_args <- function(X0, Xtau, target, S, Z, pscores, X0s, Xtaus, nj, lambda, lowlim, uplim, data_in) {
  
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
