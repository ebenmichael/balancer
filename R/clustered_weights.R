################################################################################
## Balancing weights for clustered obs studies
################################################################################



#' Re-weight groups to target population means
#' @param ind_covs n x d1 matrix of covariates for individual units
#' @param clus_covs n x d2 matrix of covariates for clusters
#' @param trt n vector of treatment assignment
#' @param clusters n vector of cluster assignments
#' @param lambda Regularization hyper parameter, default 0
#' @param icc Intraclass correlation coefficient for regularization
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
cluster_weights <- function(ind_covs, clus_covs, trt, clusters,
                            lambda = 0, icc = 0, lowlim = 0, uplim = 1,
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
    m <- length(unique(clusters[trt == 0]))
    check_data_cluster(ind_covs, clus_covs, trt, lambda, lowlim, uplim)



    # Setup the components of the QP and solve
    if(verbose) message("Creating linear term vector...")
    q <- create_q_vector_cluster(ind_covs, clus_covs, trt, m)

    if(verbose) message("Creating quadratic term matrix...")
    P <- create_P_matrix_cluster(n0, d, m, lambda / sum(trt)^2, icc)


    if(verbose) message("Creating constraint matrix...")
    constraints <- create_constraints_cluster(ind_covs, clus_covs, trt, clusters[trt == 0],
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
create_q_vector_cluster <- function(ind_covs, clus_covs, trt, m) {
    n <- nrow(ind_covs)
    n0 <- n - sum(trt)
    # concenate treated averages for each group
    q <- - c(colMeans(ind_covs[trt == 1,, drop = F]),
                     colMeans(clus_covs[trt == 1,, drop = F]))
    d <- length(q)
    q <- Matrix::sparseVector(q, (n0 + 1):(n0 + d),
                              n0 + d + m)
    return(q)
}

#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#'
#' @return P matrix
create_P_matrix_cluster <- function(n0, d, m, lambda, icc) {

  # first include a diagonal element for the auxiliary covariates X %*% gamma
  P1 <- Matrix::bdiag(Matrix::Matrix(0, n0, n0),
                      Matrix::Diagonal(d),
                      Matrix::Matrix(0, m, m))
  # Add iid variance term
  I0 <- Matrix::bdiag(Matrix::Diagonal(n0),
                      Matrix::Diagonal(d, 0),
                      Matrix::Matrix(0, m, m))
  P2 <- lambda * (1 - icc) * I0
  # add correlation within cluster term
  # cluster_mat <- Matrix::sparse.model.matrix(~ as.factor(clusters) - 1)
  P3 <- lambda * icc * Matrix::bdiag(Matrix::Diagonal(n0, 0),
                                     Matrix::Diagonal(d, 0),
                                     Matrix::Diagonal(m))
  P <- P1 + P2 + P3
  return(P)
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
create_constraints_cluster <- function(ind_covs, clus_covs, trt, clusters,
                                       lowlim, uplim, verbose) {

    d1 <- ncol(ind_covs)
    d2 <- ncol(clus_covs)
    d <- d1 + d2
    n <- nrow(ind_covs)
    n1 <- sum(trt)
    n0 <- n - n1

    cluster_mat <- Matrix::sparse.model.matrix(~ as.factor(clusters) - 1)
    m <- ncol(cluster_mat)

    if(verbose) message("\tx Sum to number of treated units constraint")
    # sum-to-n1 constraint
    A1 <- Matrix::Matrix(1, nrow = 1, ncol = n0)
    A1 <- Matrix::cbind2(A1,
                         Matrix::Matrix(0, nrow = nrow(A1), ncol = d))
    A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow = nrow(A1), ncol = m))
    l1 <- n1
    u1 <- n1
    if(verbose) message("\tx Upper and lower bounds")
    # upper and lower bounds
    A2 <- Matrix::Diagonal(n0)
    A2 <- Matrix::cbind2(A2,Matrix::Matrix(0, nrow = nrow(A2), ncol = d))
    A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = m))
    l2 <- rep(lowlim * n1, n0)
    u2 <- rep(uplim * n1, n0)


    if(verbose) message("\tx Fit weights to data")
    # constrain the auxiliary weights to be sqrt(P)'gamma
    sqrtP <- t(cbind(ind_covs, clus_covs)[trt == 0, , drop = F]) / sum(trt)
    A4 <- Matrix::cbind2(sqrtP, -Matrix::Diagonal(d))
    A4 <- Matrix::cbind2(A4, Matrix::Matrix(0, nrow = d, ncol = m))
    l4 <- rep(0, d)
    u4 <- rep(0, d)


    if(verbose) message("\tx Create cluster weights")
    # create intermediate variable that is the sum of weights in a cluster
    A5 <- Matrix::cbind2(Matrix::t(cluster_mat), Matrix::Matrix(0, nrow = m, ncol = d))
    A5 <- Matrix::cbind2(A5, -Matrix::Diagonal(m))
    l5 <- rep(0, m)
    u5 <- rep(0, m)

    if(verbose) message("\tx Combining constraints")
    A <- rbind(A1, A2, A4, A5)
    l <- c(l1, l2, l4, l5)
    u <- c(u1, u2, u4, u5)

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

#' Compute point estimate and standard error with clustered weights
#' @param y Vector of outcomes
#' @param wts Vector of weight 
#' @param trt Vector of treatment assignments
#' @param clusters Vector of cluster assignments
#' @param m1hat Vector of model predictions of E[Y(1) | covariates]
#' @param m0hat Vector of model predictions of E[Y(0) | covariates]
#' 
#' @return Data.Frame with the point estimate and standard error
#' @export
compute_cluster_se <- function(y, wts, trt, clusters, m1hat, m0hat) {

  cluster_mat <- Matrix::sparse.model.matrix(~ as.factor(clusters) - 1)
  mhat <- m1hat * trt + m0hat * (1 - trt)
  resids <- y - mhat
  n1 <- sum(trt)
  n0 <- sum(1 - trt)


  mu1 <- sum(y[trt == 1] * wts[trt == 1]) / sum(wts[trt == 1])
  mu0 <- sum(y[trt == 0] * wts[trt == 0]) / sum(wts[trt == 0])


  wtd_resids1 <- Matrix::t(cluster_mat[trt == 1, ]) %*% (y - mu1)[trt == 1]
  wtd_resids0 <- Matrix::t(cluster_mat[trt == 0, ]) %*% (wts * resids)[trt == 0]
  wtd_resids <- Matrix::t(cluster_mat) %*% (wts * resids)
  m1hat_clus <- Matrix::t(cluster_mat[trt == 1, ]) %*% (m1hat[trt == 1] - mu1)
  m0hat_clus <- Matrix::t(cluster_mat[trt == 1, ]) %*% (m0hat[trt == 1] - mu0)
  tauhat_clus <- Matrix::t(cluster_mat[trt == 1, ]) %*% ((m1hat - m0hat)[trt == 1] - (mu1 - mu0))
  se1 <- sqrt(sum(wtd_resids1^2)/sum(wts[trt == 1])^2)
  se0 <- sqrt(sum(wtd_resids0^2)/sum(wts[trt == 0])^2 +  sum(m1hat_clus^2) / sum(trt == 1)^2)
  se_tau <- sqrt(sum(wtd_resids^2)/sum(wts[trt == 1])^2 +  sum(tauhat_clus^2) / sum(trt == 1)^2)


  return(data.frame(Estimand = c("Average Outcome for Treated", "Average Counterfactual Outcome for Treated", "ATT"),
                    Estimate = c(mu1, mu0, mu1 - mu0), SE = c(se1, se0, se_tau)))
}