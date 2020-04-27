################################################################################
## Balancing weights for clustered obs studies
################################################################################



#' Re-weight groups to target population means
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param clus_covs n x d2 matrix of covariates for clusters
#' @param trt n vector of treatment assignment
#' @param lambda Regularization hyper parameter, default 0
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
#' @param verbose Whether to show messages, default T
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#'
#' @return \itemize{
#'          \item{weights }{Estimated primal weights}
#'          \item{ind_imbalance }{Imbalance in individual-level covariates covariates}
#'          \item{clus_imbalance }{Imbalance in cluster-level covariates covariates}
#'}
#' @export
cluster_weights <- function(ind_covs, clus_covs, trt, 
                            lambda = 0, lowlim = 0, uplim = 1,
                            verbose = TRUE,
                            eps_abs = 1e-5, eps_rel = 1e-5, ...) {

    # convert covariates to matrix and combine
    ind_covs <- as.matrix(ind_covs)
    clus_covs <- as.matrix(clus_covs)
    d1 <- ncol(ind_covs)
    d2 <- ncol(clus_covs)
    d <- d1 + d2
    n <- nrow(ind_covs)
    n1 <- sum(trt)
    n0 <- n - n1
    check_data_cluster(ind_covs, clus_covs, trt, lambda, lowlim, uplim)



    # Setup the components of the QP and solve
    if(verbose) message("Creating linear term vector...")
    q <- create_q_vector_cluster(ind_covs, clus_covs, trt)

    if(verbose) message("Creating quadratic term matrix...")
    P <- create_P_matrix_cluster(n0, d)

    I0 <- Matrix::bdiag(Matrix::Diagonal(n0),
                        Matrix::Diagonal(d, 0))
    P <- P + lambda * I0

    if(verbose) message("Creating constraint matrix...")
    constraints <- create_constraints_cluster(ind_covs, clus_covs, trt,
                                           lowlim, uplim, verbose)

    settings <- do.call(osqp::osqpSettings,
                        c(list(verbose = verbose,
                               eps_rel = eps_rel,
                               eps_abs = eps_abs),
                        list(...)))

    solution <- osqp::solve_osqp(P, q, constraints$A,
                                    constraints$l, constraints$u,
                                    pars = settings)
    weights <- solution$x[1:n0]
    # compute imbalance matrix
    ind_imbalance <- colMeans(ind_covs[trt == 1,, drop = F]) -
                      t(weights) %*% ind_covs[trt == 0,, drop = F] / n1
    clus_imbalance <- colMeans(clus_covs[trt == 1,, drop = F]) -
                      t(weights) %*% clus_covs[trt == 0,, drop = F] / n1
    weights <- numeric(n)
    weights[trt == 0] <- solution$x[1:n0]
    # compute overall imbalance
    return(list(weights = weights,
                ind_imbalance = ind_imbalance,
                clus_imbalance = clus_imbalance))
}


#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param clus_covs n x d2 matrix of covariates for clusters
#' @param trt n vector of treatment assignment
#' 
#' @return q vector
create_q_vector_cluster <- function(ind_covs, clus_covs, trt) {
    n <- nrow(ind_covs)
    n0 <- n - sum(trt)
    # concenate treated averages for each group
    q <- - c(colSums(ind_covs[trt == 1,, drop = F]),
                     colSums(clus_covs[trt == 1,, drop = F]))
    d <- length(q)
    q <- Matrix::sparseVector(q, (n0 + 1):(n0 + d),
                              n0 + d)
    return(q)
}

#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#'
#' @return P matrix
create_P_matrix_cluster <- function(n0, d) {
    return(Matrix::bdiag(Matrix::Matrix(0, n0, n0),
                         Matrix::Diagonal(d)))
}



#' Create the constraints for QP: l <= Ax <= u
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param clus_covs n x d2 matrix of covariates for clusters
#' @param trt n vector of treatment assignment
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#' @param verbose Whether to show messages
#'
#' @return A, l, and u
create_constraints_cluster <- function(ind_covs, clus_covs, trt,
                                       lowlim, uplim, verbose) {

    d1 <- ncol(ind_covs)
    d2 <- ncol(clus_covs)
    d <- d1 + d2
    n <- nrow(ind_covs)
    n1 <- sum(trt)
    n0 <- n - n1

    if(verbose) message("\tx Sum to number of treated units constraint")
    # sum-to-n1 constraint for each group
    A1 <- Matrix::Matrix(1, nrow = 1, ncol = n0)
    A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow = nrow(A1), ncol = d))
    l1 <- n1
    u1 <- n1
    if(verbose) message("\tx Upper and lower bounds")
    # upper and lower bounds
    A2 <- Matrix::Diagonal(n0)
    A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = d))
    l2 <- rep(lowlim * n1, n0)
    u2 <- rep(uplim * n1, n0)


    if(verbose) message("\tx Fit weights to data")
    # constrain the auxiliary weights to be sqrt(P)'gamma
    sqrtP <- t(cbind(ind_covs, clus_covs)[trt == 0, , drop = F])
    A4 <- Matrix::cbind2(sqrtP, -Matrix::Diagonal(d))
    l4 <- rep(0, d)
    u4 <- rep(0, d)

    if(verbose) message("\tx Combining constraints")
    A <- rbind(A1, A2, A4)
    l <- c(l1, l2, l4)
    u <- c(u1, u2, u4)

    return(list(A = A, l = l, u = u))
}



#' Check that data is in right shape and hyparameters are feasible
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param clus_covs n x d2 matrix of covariates for clusters
#' @param trt n vector of treatment assignment
#' @param lambda Regularization hyper parameter
#' @param lowlim Lower limit on weights, default 0
#' @param uplim Upper limit on weights, default 1
check_data_cluster <- function(ind_covs, clus_covs, trt, lambda, lowlim, uplim) {


    # NA checks
    if(any(is.na(ind_covs))) {
        stop("Individuals covariate matrix contains NA values.")
    }

    if(any(is.na(clus_covs))) {
        stop("Clusters covariate matrix contains NA values.")
    }

    if(any(! trt %in% c(0,1))) {
        stop("Treatment must be (0,1)")
    }

    #dimension checks
    n1 <- nrow(ind_covs)
    n2 <- nrow(clus_covs)
    d1 <- ncol(ind_covs)
    d2 <- ncol(clus_covs)

    if(n1 != n2) {
      stop("Individuals covariate matrix and clusters covariate matrix do not have the same number of rows")
    }

    # hyerparameters are feasible
    if(lambda < 0) {
        stop("lambda must be >= 0")
    }
    if(lowlim > uplim) {
        stop("Lower threshold must be lower than upper threshold")
    }
}