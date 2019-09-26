################################################################################
## Wrapper to standardize to target means
################################################################################

#' Re-weight groups to target population means
#' @param X n x d matrix of covariates
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators with J levels
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_data Whether to return the objective matrix and vector and constraints, default T
#'
#' @return \itemize{
#'          \item{weights }{Estimated primal weights as an n x J matrix}
#'          \item{imbalance }{Imbalance in covariates as a d X J matrix}
#'          \item{data_out }{List containing elements of QP min 0.5 x'Px + q'x st l <= Ax <= u \itemize{
#'                  \item{P, q}{}
#'                  \item{constraints }{A, l , u}
#'}}}
#' @export
standardize <- function(X, target, Z, lambda, lowlim = 0, uplim = 1, 
                        data_in = NULL, verbose = TRUE, return_data = TRUE) {

    Z_factor <- as.factor(Z)
    unique_Z <- levels(Z_factor)
    J <- length(unique_Z)
    # dimension of auxiliary weights
    aux_dim <- J * ncol(X)
    n <- nrow(X)

    # split matrix by targets
    Xz <- split.data.frame(X, Z_factor)
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
        P <- data_in$P + lambda * Matrix::Diagonal(nrow(X))
    }
    I0 <- Matrix::bdiag(Matrix::Diagonal(nrow(X)), Matrix::Diagonal(aux_dim, 0))
    P <- P + lambda * I0

    if(verbose) message("Creating constraint matrix...")
    if(is.null(data_in$constraints)) {
        constraints <- create_constraints(Xz, target, Z, lowlim, uplim, verbose)
    } else {
        constraints <- data_in$constraints
    }


    if(verbose) {
        settings <- osqp::osqpSettings(verbose = TRUE, eps_abs = 1e-6,
                                   eps_rel = 1e-6, max_iter = 5000)
    } else {
        settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-6,
                                   eps_rel = 1e-6, max_iter = 5000)
    }
    solution <- osqp::solve_osqp(P, q, constraints$A,
                                 constraints$l, constraints$u, pars = settings)


    # convert weights into a matrix

    nj <- sapply(1:J, function(j) nrow(Xz[[j]]))
    weights <- Matrix::Matrix(0, ncol = J, nrow = n)

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

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param Xz list of J n x d matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param aux_dim Dimension of auxiliary weights
#'
#' @return q vector
create_q_vector <- function(Xz, target, aux_dim) {
    q <- c(do.call(rbind, Xz) %*% target)
    q <- Matrix::sparseVector(q, 1:length(q),
                              length(q) + aux_dim)
    return(q)
}


#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#' @param X n x d matrix of covariates
#' @param Z Vector of group indicators
#'
#' @return P matrix
create_P_matrix <- function(n, aux_dim) {
    return(Matrix::Diagonal(n + aux_dim))
}

#' Create the constraints for QP: l <= Ax <= u
#' @param Xz list of J n x d matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#'
#' @return A, l, and u
create_constraints <- function(Xz, target, Z, lowlim, uplim, verbose) {

    J <- length(Xz)
    nj <- sapply(1:J, function(j) nrow(Xz[[j]]))
    d <- ncol(Xz[[1]])
    cumsum_nj <- cumsum(c(1, nj))
    n <- sum(nj)
    Xzt <- lapply(Xz, t)

    # dimension of auxiliary weights
    aux_dim <- J * d

    if(verbose) message("\tx Sum to one constraint")
    # sum-to-one constraint for each group
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

    if(verbose) message("\tx Mantain overall population mean")
    # Constrain the overall mean to be equal to the target
    A3 <- do.call(cbind, lapply(Xzt, function(x) x * ncol(x)))

    A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))

    l3 <- n * target
    u3 <- n * target

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


standardize_old <- function(X, target,
                        Z=NULL, 
                        link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                        regularizer=c("ridge", "l1", "grpl1", "l2", "linf", "nuc",
                                      "l1_all", "l1_nuc"),
                        lambda=NULL, nlambda=20, lambda.min.ratio=1e-3,
                        normalized=TRUE, alpha=1,
                        opts=list()) {
    #' Re-weight groups to target population means
    #' @param X n x d matrix of covariates
    #' @param target Vector of population moments to weight towards
    #' @param Z Vector of site indicators
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param lambda Regularization hyperparameter
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param interact Whether to interact group and individual level covariates
    #' @param normalized Whether to normalize the weights
    #' @param alpha Elastic net parameter \eqn{\frac{1-\alpha}{2}\|\beta\|_2^2 + \alpha\|\beta\|_1}, defaults to 1
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}
    #'          \item{alpha }{Elastic net parameter}}
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}
    #' @export


    ## map string args to actual params
    suppressWarnings(params <- map_to_param(X, link, regularizer, NULL, normalized, Q, alpha))
    weightfunc <- params[[1]]
    weightptr <- params[[2]]
    proxfunc <- params[[3]]
    balancefunc <- params[[4]]
    prox_opts <- params[[5]]


    ipw_weights = matrix(1/nrow(X), nrow(X), 1)
    
    ## add intercept
    if(normalized) {
        X <- cbind(1, X)
        target <- c(1, target)
    }

    out <- standardize_(X, target, Z, weightfunc, weightptr, proxfunc,
                        balancefunc, lambda, nlambda, lambda.min.ratio,
                        ipw_weights, opts, prox_opts)
    return(out)
}

standardize_  <- function(X, target, Z, weightfunc, weightfunc_ptr,
                            proxfunc, balancefunc, lambda, nlambda,
                            lambda.min.ratio, ipw_weights, 
                          opts=list(), prox_opts) {
    #' Internal function to reweight to target population means
    #' @param X n x d matrix of covariates
    #' @param target Vector of population means to re-weight to
    #' @param Z Vector of hierarchical factor indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param balancefunc Balance criterion measure
    #' @param lambda Regularization hyper parameter
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}}
    #' @param prox_opts List of additional arguments for prox    
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}
    
    ## if no subgroups, put everything into one group
    if(is.null(Z)) {
        Z <- rep(0, nrow(X))
    }


    
    ## get the distinct group labels
    grps <- sort(unique(Z))
    m <- length(grps)

    n <- dim(X)[1]
    d <- dim(X)[2]
    
    target <- as.matrix(sapply(grps,
                  function(k) target))


    loss_opts = list(X=X,
                     target=target,
                     z=Z,
                     weight_func=weightfunc_ptr,
                     ipw_weights=ipw_weights,
                     n_groups=m,
                     dim=d
                     )


    ## initialize at 0
    init = matrix(0, nrow=d, ncol=(m))
    
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))

    ## if hyperparam is NULL, start from reference weights and decrease
    if(is.null(lambda)) {
        lam0 <- balancefunc(balancing_grad_multilevel(init, loss_opts))
        lam1 <- lam0 * lambda.min.ratio
        ## decrease on log scale
        lambda <- exp(seq(log(lam0), log(lam1), length.out=nlambda))
    }


    ## collect results
    out <- list()
    out$theta <- matrix(,nrow=d, ncol=length(lambda))
    out$imbalance <- matrix(,nrow=d, ncol=length(lambda))    
    out$weights <- matrix(0, nrow=n, ncol=length(lambda))
    out$weightfunc <- weightfunc


    ## with multiple hyperparameters do warm starts        
    prox_opts <- c(prox_opts,
                  list(lam=1))


    apgout <- apg_warmstart(make_balancing_grad_standardize(),
                            proxfunc, loss_opts, prox_opts,
                            lambda,
                            opts$x, opts$max_it, opts$eps,
                            opts$alpha, opts$beta, opts$accel, opts$verbose)

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    weights_from_theta <- function(th) {
        weights <- matrix(0, nrow=n, ncol=m)
        for(i in 1:m) {
            k = grps[i]
            weights[(Z == k), i] <-
                weightfunc(X[Z == k, , drop=FALSE],
                           th[,i,drop=FALSE],
                           ipw_weights[Z == k,,drop=FALSE])
        }
        weights
    }
    
    weights <- lapply(out$theta, weights_from_theta)

    out$weights <- if(length(weights)==1) weights[[1]] else weights

    
    ## The final imbalance    
    imbalance <- lapply(theta, function(th) balancing_grad_standardize(th, loss_opts))
    out$imbalance <- if(length(imbalance) == 1) imbalance[[1]] else imbalance

    out$theta <- if(length(out$theta) == 1) out$theta[[1]] else out$theta

    out$lambda <- lambda
    return(out)
}