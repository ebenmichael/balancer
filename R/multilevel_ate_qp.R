#' Re-weight treated and control sub-groups to sub-group means
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment assignments
#' @param Z Vector of group indicators with J levels
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param exact_global Whether to enforce exact balance for overall population
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
multilevel_ate_qp <- function(X, trt, Z, lambda = 0, lowlim = 0, uplim = 1,
                        scale_sample_size = T, exact_global = T,
                        verbose = TRUE,
                        eps_abs = 1e-5, eps_rel = 1e-5, ...) {
    # # convert X to a matrix
    # X <- as.matrix(X)

    # # split data by factor and compute average within group
    # Z_factor <- as.factor(Z)
    # targetz <- lapply(split.data.frame(X, Z_factor), colMeans)
    # target_propz <- sapply(split(trt, Z_factor), length) / nrow(X)
    
    # # split into treatment and control
    # Xz0 <- split.data.frame(X[trt == 0, ], Z_factor)
    # Xz1 <- split.data.frame(X[trt == 1, ], Z_factor)

    # check_data_multi(X, trt, Z, split.data.frame(X, Z_factor), lambda, lowlim, uplim)


    # weights_ctrl <- l2_balance_internal(Xz0, targetz, lambda, lowlim, uplim,
    #                                     scale_sample_size, exact_global,
    #                                     target_propz, verbose,
    #                                     eps_abs, eps_rel, ...)
    n <- nrow(X)
    all_trt <- cbind(numeric(n), rep(1, n))
    colnames(all_trt) <- c(0, 1)
    all_ctrl <- cbind(rep(1, n), numeric(n))
    colnames(all_ctrl) <- c(0, 1)


    # control weights
    ctrl_wts <- stochastic_int(X, trt, all_ctrl,
                               Z, lambda, lowlim, uplim,
                               exact_global, verbose,
                               eps_abs, eps_rel, ...)
    
    # treated weights
    # control weights
    trt_wts <- stochastic_int(X, trt, all_trt,
                               Z, lambda, lowlim, uplim,
                               exact_global, verbose,
                               eps_abs, eps_rel, ...)
    
    # combine into one set of weights

    full_weights <- numeric(n)
    full_weights[trt == 1] <- trt_wts$weights[trt == 1]
    full_weights[trt == 0] <- ctrl_wts$weights[trt == 0]

    return(list(
      weights = full_weights,
      imbalance_ctrl = ctrl_wts$imbalance,
      global_imbalance_ctrl = ctrl_wts$global_imbalance,
      imbalance_trt = trt_wts$imbalance,
      global_imbalance_trt = trt_wts$global_imbalance
    ))

}