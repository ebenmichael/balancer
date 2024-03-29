################################################################################
## Multilevel balancing weights
################################################################################

#' Re-weight control sub-groups to treated sub-group means
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
multilevel_qp <- function(X, trt, Z, lambda = 0, lowlim = 0, uplim = 1,
                        scale_sample_size = T, exact_global = T,
                        verbose = TRUE,
                        eps_abs = 1e-5, eps_rel = 1e-5, ...) {

    # convert X to a matrix
    X <- as.matrix(X)

    # split data and treatment by factor
    Z_factor <- as.factor(Z)
    Xz <- split.data.frame(X, Z_factor)
    trtz <- split(trt, Z)

    check_data_multi(X, trt, Z, Xz, lambda, lowlim, uplim)


    unique_Z <- levels(Z_factor)
    J <- length(unique_Z)
    # dimension of auxiliary weights
    aux_dim <- J * ncol(X)
    n <- nrow(X)




    idxs <- split(1:nrow(X), Z_factor)



    # Setup the components of the QP and solve
    if(verbose) message("Creating linear term vector...")
    q <- create_q_vector_multi(Xz, trtz)

    if(verbose) message("Creating quadratic term matrix...")
    P <- create_P_matrix_multi(n, aux_dim)

    I0 <- create_I0_matrix_multi(Xz, scale_sample_size, n, aux_dim)
    P <- P + lambda * I0

    if(verbose) message("Creating constraint matrix...")
    constraints <- create_constraints_multi(Xz, trtz, lowlim, 
                                            uplim, exact_global, verbose)

    settings <- do.call(osqp::osqpSettings,
                        c(list(verbose = verbose,
                               eps_rel = eps_rel,
                               eps_abs = eps_abs),
                        list(...)))

    solution <- osqp::solve_osqp(P, q, constraints$A,
                                    constraints$l, constraints$u,
                                    pars = settings)

    # convert weights into a matrix
    nj <- sapply(1:J, function(j) nrow(Xz[[j]]))
    weights <- numeric(n)

    if(verbose) message("Reordering weights...")
    cumsumnj <- cumsum(c(1, nj))
    for(j in 1:J) {
        weights[idxs[[j]]] <- solution$x[cumsumnj[j]:(cumsumnj[j + 1] - 1)]
    }

    # compute imbalance matrix
    n1j <- sapply(trtz, sum)
    imbalance <- vapply(1:J,
                        function(j) {
                            target <- colMeans(Xz[[j]][trtz[[j]] == 1, , drop = F])
                            target - t(X[idxs[[j]], ]) %*% weights[idxs[[j]]] / n1j[[j]]
                        },
                        numeric(ncol(X)))

    # compute overall imbalance
    global_imbal <- colSums(t(imbalance) * n1j) / sum(n1j)
    global_imbal <-  colMeans(X[trt == 1,, drop = F]) - t(X) %*% weights / sum(n1j)
    return(list(weights = weights,
                imbalance = imbalance,
                global_imbalance = global_imbal))

}

#' Create diagonal regularization matrix
#' @param Xz list of J n x d matrices of covariates split by group
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param n Total number of units
#' @param aux_dim Dimension of auxiliary weights
create_I0_matrix_multi <- function(Xz, scale_sample_size, n, aux_dim) {

    if(scale_sample_size) {
        # diagonal matrix 1 / n_j for each group j
        subdiags <- lapply(Xz,
                        function(x) Matrix::Diagonal(nrow(x), 1 / nrow(x)))
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
create_q_vector_multi <- function(Xz, trtz) {
    J <- length(Xz)
    d <- ncol(Xz[[1]])
    n <- Reduce(`+`, lapply(Xz, nrow))
    aux_dim <- J * d
    # concenate treated averages for each group
    q <- - do.call(c,
                 lapply(1:J,
                    function(j) colSums(Xz[[j]][trtz[[j]] == 1, ,drop = F])
                       )
                )

    q <- Matrix::sparseVector(q, (n + 1):(n + aux_dim),
                              n + aux_dim)
    return(q)
}


#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#' @param X n x d matrix of covariates
#' @param Z Vector of group indicators
#'
#' @return P matrix
create_P_matrix_multi <- function(n, aux_dim) {
    return(Matrix::bdiag(Matrix::Matrix(0, n, n),
                         Matrix::Diagonal(aux_dim)))
}

# #' Get a set of uniform weights for initialization
# #' @param Xz list of J n x d matrices of covariates split by group
# #' 
# get_uniform_weights <- function(Xz) {

#     # uniform weights for each group
#     uniw <- do.call(c, lapply(Xz, function(x) rep(1 / nrow(x), nrow(x))))

#     # transformed auxiliary uniform weights
#     sqrtP <- Matrix::bdiag(lapply(Xz, t))
#     aux_uniw <- as.numeric(sqrtP %*% uniw)
#     return(c(uniw, aux_uniw))
# }

#' Create the constraints for QP: l <= Ax <= u
#' @param Xz list of J n x d matrices of covariates split by group
#' @param trtz list of treatment assignment vectors
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param exact_global Boolean whether to include an exact global constraint
#' @param verbose Boolean whether to display progress
#'
#' @return A, l, and u
create_constraints_multi <- function(Xz, trtz, lowlim, uplim, exact_global, verbose) {

    J <- length(Xz)
    n0j <- sapply(1:J, function(j) nrow(Xz[[j]]))
    n1j <- sapply(trtz, sum)

    d <- ncol(Xz[[1]])
    n <- Reduce(`+`, lapply(Xz, nrow))
    Xzt <- lapply(Xz, t)

    # dimension of auxiliary weights
    aux_dim <- J * d

    if(verbose) message("\tx Sum to one constraint")
    # sum-to-n1j constraint for each group
    A1 <- Matrix::t(Matrix::bdiag(lapply(Xz, function(x) rep(1, nrow(x)))))
    A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = aux_dim))
    l1 <- n1j
    u1 <- n1j
    if(verbose) message("\tx Upper and lower bounds")
    # upper and lower bounds
    A2 <- Matrix::Diagonal(n)
    A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
    l2 <- Reduce(c, sapply(1:J, function(j) rep(lowlim * n1j[[j]], n0j[j])))
    u2 <- Reduce(c, sapply(1:J, function(j) rep(uplim * n1j[[j]], n0j[j])))

    if(exact_global) {
        if(verbose) message("\tx Enforce exact global balance")
        # Constrain the overall mean to be equal to the target
        A3 <- do.call(cbind, lapply(1:J, function(j) Xzt[[j]]))
        A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
        trt_sum <- Reduce(`+`,
            lapply(1:J,
                function(j) colSums(Xz[[j]][trtz[[j]] == 1, , drop = F])))
        l3 <- trt_sum
        u3 <- trt_sum
        
    } else {
        if(verbose) message("\t(SKIPPING) Enforce exact global balance")
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

    if(verbose) message("\tx Constrain treated weights to be zero")
    # zero out treated units
    A5 <- Matrix::bdiag(lapply(trtz, function(x) Matrix::Diagonal(x = x)))
    A5 <- Matrix::cbind2(A5, Matrix::Matrix(0, nrow = nrow(A5), ncol = aux_dim))
    l5 <- numeric(n)
    u5 <- numeric(n)

    if(verbose) message("\tx Combining constraints")
    A <- rbind(A1, A2, A3, A4, A5)
    l <- c(l1, l2, l3, l4, l5)
    u <- c(u1, u2, u3, u4, u5)

    return(list(A = A, l = l, u = u))
}



#' Check that data is in right shape and hyparameters are feasible
#' @param X n x d matrix of covariates
#' @param Z Vector of group indicators with J levels
#' @param Xz list of J n x d matrices of covariates split by group
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
check_data_multi <- function(X, trt, Z, Xz, lambda, lowlim, uplim) {


    # NA checks
    if(any(is.na(X))) {
        stop("Covariate matrix X contains NA values.")
    }

    if(any(! trt %in% c(0,1))) {
        stop("Treatment must be (0,1)")
    }

    if(any(is.na(Z))) {
        stop("Grouping vector Z contains NA values.")
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


    # hyerparameters are feasible
    if(lambda < 0) {
        stop("lambda must be >= 0")
    }
    if(lowlim > uplim) {
        stop("Lower threshold must be lower than upper threshold")
    }
    if(lowlim > 1 / max(nj)) {
        stop("Lower threshold must be lower than 1 / size of largest group")
    }
    if(uplim < 1 / min(nj)) {
        stop("Upper threshold must be higher than 1 / size of smallest group")
    }

}