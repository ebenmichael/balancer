################################################################################
## Wrapper to standardize to target means
################################################################################

#' Re-weight groups to target population means
#' @param X n x d matrix of covariates
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators with J levels
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_data Whether to return the objective matrix and vector and constraints, default T
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
standardize <- function(X, target, Z, lambda = 0, lowlim = 0, uplim = 1,
                        scale_sample_size = T,
                        data_in = NULL, verbose = TRUE, return_data = TRUE,
                        exact_global = T, init_uniform = F,
                        eps_abs = 1e-5, eps_rel = 1e-5, ...) {

    # convert X to a matrix
    X <- as.matrix(X)

    # split matrix by targets
    Z_factor <- as.factor(Z)
    Xz <- split.data.frame(X, Z_factor)

    # ensure that target is a vector
    target <- c(target)

    check_data(X, target, Z, Xz, lambda, lowlim, uplim, data_in)


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
        constraints <- create_constraints(Xz, target, Z, lowlim,
                                          uplim, exact_global, verbose)
    } else {
        constraints <- data_in$constraints
        constraints$l[(J + 1): (J + n)] <- lowlim
        constraints$u[(J + 1): (J + n)] <- uplim
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


#' Create diagonal regularization matrix
#' @param Xz list of J n x d matrices of covariates split by group
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param n Total number of units
#' @param aux_dim Dimension of auxiliary weights
create_I0_matrix <- function(Xz, scale_sample_size, n, aux_dim) {

    if(scale_sample_size) {
        # diagonal matrix n_j / n for each group j
        subdiags <- lapply(Xz,
                        function(x) Matrix::Diagonal(nrow(x), nrow(x)))
        I0 <- Matrix::bdiag(subdiags)
    } else {
        # all diagonal entries are 1
        I0 <- Matrix::Diagonal(n)
    }
    I0 <- Matrix::bdiag(I0, Matrix::Diagonal(aux_dim, 0))
    return(I0)
}

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param Xz list of J n x d matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param aux_dim Dimension of auxiliary weights
#'
#' @return q vector
create_q_vector <- function(Xz, target, aux_dim) {
    q <- -c(do.call(rbind, Xz) %*% target)
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
    return(Matrix::bdiag(Matrix::Diagonal(n, 0), Matrix::Diagonal(aux_dim, 1)))
}

#' Get a set of uniform weights for initialization
#' @param Xz list of J n x d matrices of covariates split by group
#'
get_uniform_weights <- function(Xz) {

    # uniform weights for each group
    uniw <- do.call(c, lapply(Xz, function(x) rep(1 / nrow(x), nrow(x))))

    # transformed auxiliary uniform weights
    sqrtP <- Matrix::bdiag(lapply(Xz, t))
    aux_uniw <- as.numeric(sqrtP %*% uniw)
    return(c(uniw, aux_uniw))
}

#' Create the constraints for QP: l <= Ax <= u
#' @param Xz list of J n x d matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#'
#' @return A, l, and u
create_constraints <- function(Xz, target, Z, lowlim, uplim, exact_global, verbose) {

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

    if(exact_global) {
        if(verbose) message("\tx Mantain overall population mean")
        # Constrain the overall mean to be equal to the target
        A3 <- do.call(cbind, lapply(Xzt, function(x) x * ncol(x)))

        A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))

        l3 <- n * target
        u3 <- n * target
    } else {
        if(verbose) message("\t(SKIPPING) Mantain overall population mean")
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



#' Check that data is in right shape and hyparameters are feasible
#' @param X n x d matrix of covariates
#' @param target Vector of population means to re-weight to
#' @param Z Vector of group indicators with J levels
#' @param Xz list of J n x d matrices of covariates split by group
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param data_in Optional list containing pre-computed objective matrix/vector and constraints (without regularization term)
#' @param verbose Whether to show messages, default T
#' @param return_data Whether to return the objective matrix and vector and constraints, default T
check_data <- function(X, target, Z, Xz, lambda, lowlim, uplim, data_in) {


    # NA checks
    if(any(is.na(X))) {
        stop("Covariate matrix X contains NA values.")
    }

    if(any(is.na(Z))) {
        stop("Grouping vector Z contains NA values.")
    }

    if(any(is.na(target))) {
        stop("Target vector contains NA values.")
    }

    #dimension checks
    n <- nrow(X)
    d <- ncol(X)
    J <- length(Xz)
    aux_dim <- d * J
    nj <- as.numeric(lapply(Xz, nrow))

    if(length(Z) != n) {
        stop("The number of rows in covariate matrix X (", n,
             ") does not equal the dimension of and grouping vector Z (",
             length(Z), ").")
    }

    if(sum(nj) != n) {
        stop("Implied number of weights (", sum(nj),
             ") does not equal number of units (", n, ").")
    }

    if(length(target) != d) {
        stop("Target dimension (", length(target),
             ") is not equal to data dimension (", d, ").")
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
        if(length(data_in$constraints$l) != J + n + d +aux_dim) {
            stop("data_in$constraints$l must have dimension ",
                 J + d + n + aux_dim)
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
