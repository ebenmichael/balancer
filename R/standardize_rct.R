################################################################################
## Wrapper to standardize multisite RCTs to a target level
################################################################################

#' Re-weight groups to target population means
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment assignments
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators with J levels
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_data Whether to return the objective matrix and vector and constraints, default T
#' @param constrain_trt_prop Whether to constrain the treated proportions within each site to be equal to the observed proportion
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
standardize_rct <- function(X, trt, target, Z,
                            lambda = 0, lowlim = 0, uplim = 1,
                            scale_sample_size = T,
                            data_in = NULL, verbose = TRUE, return_data = TRUE,
                            constrain_trt_prop = TRUE,
                            exact_global = F, init_uniform = F,
                            eps_abs = 1e-5, eps_rel = 1e-5, ...) {

    # convert X to a matrix
    X <- as.matrix(X)

    # split matrix by targets
    Z_factor <- as.factor(Z)
    Xz <- split.data.frame(X, Z_factor)

    # ensure that target is a vector
    target <- c(target)

    check_data_rct(X, trt, target, Z, Xz, lambda, lowlim, uplim, data_in)


    unique_Z <- levels(Z_factor)
    J <- length(unique_Z)
    # dimension of auxiliary weights
    aux_dim <- J * ncol(X)
    n <- nrow(X)




    idxs <- split(1:nrow(X), Z_factor)



    # Setup the components of the QP and solve
    if(verbose) message("Creating linear term vector...")
    if(is.null(data_in$q)) {
        q <- create_q_vector(Xz, target, aux_dim)
    } else {
        q <- data_in$q
    }


    if(verbose) message("Creating quadratic term matrix...")
    if(is.null(data_in$P)) {
        P <- create_P_matrix(n, aux_dim)
    } else {
        P <- data_in$P
    }
    I0 <- create_I0_matrix(Xz, scale_sample_size, n, aux_dim)
    P <- P + lambda * I0

    if(verbose) message("Creating constraint matrix...")
    if(is.null(data_in$constraints)) {
        constraints <- create_constraints_rct(Xz, trt, target, Z, lowlim,
                                              uplim, exact_global,
                                              constrain_trt_prop, verbose)
    } else {
        constraints <- data_in$constraints
        constraints$l[(2 * J + 1): (2 * J + n)] <- lowlim
        constraints$u[(2 * J + 1): (2 * J + n)] <- uplim
    }


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

    if(init_uniform) {
        if(verbose) message("Initializing with uniform weights")
        # initialize with uniform weights
        unifw <- get_uniform_weights(Xz)
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

    nj <- sapply(1:J, function(j) nrow(Xz[[j]]))
    weights <- matrix(0, ncol = J, nrow = n)

    if(verbose) message("Reordering weights...")
    cumsumnj <- cumsum(c(1, nj))
    for(j in 1:J) {
        weights[idxs[[j]], j] <- solution$x[cumsumnj[j]:(cumsumnj[j + 1] - 1)]
    }

    # compute imbalance matrix
    imbalance <- as.matrix(target - t(X) %*% weights)

    if(return_data) {
        data_out <- list(P = P  - lambda * I0,
                         q = q, constraints = constraints)
    } else {
        data_out <- NULL
    }

    return(list(weights = weights, imbalance = imbalance, data_out = data_out))

}


#' Create the constraints for QP: l <= Ax <= u
#' @param Xz list of J n x d matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#'
#' @return A, l, and u
create_constraints_rct <- function(Xz, trt, target, Z, 
                                   lowlim, uplim, exact_global,
                                   constrain_trt_prop, verbose) {
    J <- length(Xz)
    d <- ncol(Xz[[1]])
    # dimension of auxiliary weights
    aux_dim <- J * d

    # compute constraints then add on treatment proportion constraint
    consts <- create_constraints(Xz, target, Z, lowlim, uplim,
                                 exact_global, verbose)


    if(constrain_trt_prop) {
        # split treatment vector
        trtz <- split(trt, Z)
        if(verbose) message("\tx Keep treatment proportion unchanged")
        A0 <- Matrix::t(Matrix::bdiag(trtz))
        A0 <- Matrix::cbind2(A0, 
                            Matrix::Matrix(0, nrow=nrow(A0), ncol = aux_dim))
        # treated proportions in each group
        p1z <- sapply(trtz, mean)
        l0 <- p1z
        u0 <- p1z
        consts$A <- rbind(A0, consts$A)
        consts$l <- c(l0, consts$l)
        consts$u <- c(u0, consts$u)
    }


    return(consts)
}


#' Check that data is in right shape and hyparameters are feasible
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment assignments
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators with J levels
#' @param Xz list of J n x d matrices of covariates split by group
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_data Whether to return the objective matrix and vector and constraints, default T
check_data_rct <- function(X, trt, target, Z, Xz, lambda, 
                           lowlim, uplim, data_in) {

    if(any(is.na(trt))) {
        stop("Treatment vector contains NA values.")
    }

    if(length(trt) != nrow(X)) {
        stop("The number of rows in covariate matrix (", n,
             ") does not equal the dimension of and treatment vector (",
             length(trt), ").")
    }
    
    check_data(X, target, Z, Xz, lambda, lowlim, uplim, data_in) 
}