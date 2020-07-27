################################################################################
## Wrapper to standardize to target means
################################################################################

#' Re-weight groups to target population means
#' @param X0 n x d0 matrix of untransformed covariates defining the mean control response function
#' @param Xtau n x dtau matrix of transformed covariates defining the mean treatment effect function
#' @param target Vector of population means to re-weight to
#' @param S Numeric vector of site indicators with J levels
#' @param Z Numeric vector of treatment indicators with 2 levels
#' @param pscores Numeric vector of propensity scores
#' @param kernel0 Kernel for control outcome covariates, default is the inner product
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
#' @param gc boolean indicating whether to garbage collect between operations
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
standardize_treatment_hybrid <- function(X0, Xtau, target, S, Z, pscores,
                                         kernel0 = kernlab::vanilladot(), lambda = 0,
                                         lowlim = 0, uplim = nrow(X0), scale_sample_size = F,
                                         data_in = NULL, verbose = TRUE,
                                         return_program = TRUE,
                                         init_uniform = F, eps_abs = 1e-5,
                                         eps_rel = 1e-5, gc, ...) {

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
  check_args_treatment_kernel(X0, Xtau, target, S, Z, pscores, X0s, Xtaus,
                              as.numeric(nj), lambda, lowlim, uplim, data_in)

  # create propensity score multipliers and split by site
  pro_trt <- Z / (pscores * nj[S_factor])
  pro_trt_split <- split(pro_trt, S_factor)
  pro_ctr <- (1 - Z) / ((1 - pscores) * nj[S_factor])
  pro_ctr_split <- split(pro_ctr, S_factor)

  # construct linear term vector
  if(verbose) message("Creating linear term vector...")
  tic("q vector")
  if(is.null(data_in$q)) {
    q <- create_q_vector_treatment(Xtaus, pro_trt_split, target, 0)
  } else {
    q <- data_in$q
  }
  toc(log = TRUE)

  # construct quadratic term matrix
  if(verbose) message("Creating quadratic term matrix...")
  if(is.null(data_in$P)) {
    P <- create_P_matrix_treatment_kernel(n, X0s, Xtaus, kernel0, kernlab::vanilladot(),
                                          pro_trt, pro_ctr, S_factor, gc)
  } else {
    P <- data_in$P
  }
  tic("create IO matrix)")
  I0 <- create_I0_matrix_treatment(pro_trt_split, pro_ctr_split,
                                   scale_sample_size, nj, n, 0)
  toc(log = TRUE)
  tic("P + I0")
  P <- P + lambda * I0
  toc(log = TRUE)

  # construct constraint matrix
  if(verbose) message("Creating constraint matrix...")
  tic("constraint matrix creation")
  if(is.null(data_in$constraints)) {
    constraints <- create_constraints_treatment_kernel(X0s, Xtaus, Z, S_factor,
                                                       pro_trt_split, pro_ctr_split,
                                                       lowlim, uplim, verbose)
  } else {
    constraints <- data_in$constraints
    constraints$l[(J + 1):(J + n)] <- lowlim
    constraints$u[(J + 1):(J + n)] <- uplim
  }
  toc(log = TRUE)

  # set optimization settings
  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                             eps_rel = eps_rel,
                             eps_abs = eps_abs),
                        list(...)))

  # solve optimization problem (possibly with uniform weights)
  tic("solve QP")
  if(init_uniform) {
    if(verbose) message("Initializing with uniform weights")
    unifw <- get_uniform_weights_treatment_kernel(nj)
    obj <- osqp::osqp(P, q, constraints$A,
                      constraints$l, constraints$u, pars = settings)
    obj$WarmStart(x = unifw)
    solution <- obj$Solve()
  } else {
    solution <- osqp::solve_osqp(P, q, constraints$A,
                                 constraints$l, constraints$u,
                                 pars = settings)
  }
  toc(log = TRUE)

  tic("post-processing")
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
  toc(log = TRUE)

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
