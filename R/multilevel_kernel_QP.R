################################################################################
## Multilevel balancing weights
################################################################################

#' Re-weight control sub-groups to treated sub-group with kernel imbalance
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment assignments
#' @param Z Vector of group indicators with J levels
#' @param kernel Kernel to compute balance measure with, default is the inner product
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
multilevel_kernel_qp <- function(X, trt, Z,
                                kernel = kernlab::vanilladot(),
                                lambda = 0, lowlim = 0, uplim = 1,
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

    kern_mats <- compute_kernel(Xz, kernel)


    # Setup the components of the QP and solve
    if(verbose) message("Creating linear term vector...")
    q <- create_q_vector_multi_kern(kern_mats, trtz)


    if(verbose) message("Creating quadratic term matrix...")
    P <- create_P_matrix_multi_kern(kern_mats)

    I <- create_scaled_I_matrix(Xz, scale_sample_size)
    P <- P + lambda * I

    if(verbose) message("Creating constraint matrix...")
    constraints <- create_constraints_multi_kern(Xz, trtz, lowlim, 
                                            uplim, exact_global, 
                                            verbose)

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




#' compute kernel matrix
compute_kernel <- function(Xz, kernel) {

    # individual kernel matrices 
    kern_mats <- lapply(Xz,
                        function(x) kernlab::kernelMatrix(kernel, x))
    return(kern_mats)
    
}



#' Create diagonal regularization matrix
#' @param Xz list of J n x d matrices of covariates split by group
#' @param scale_sample_size Whether to scale the dispersion penalty by the sample size of each group, default T
#' @param n Total number of units
#' @param aux_dim Dimension of auxiliary weights
create_scaled_I_matrix <- function(Xz, scale_sample_size) {

    if(scale_sample_size) {
        # diagonal matrix n_j / n for each group j
        subdiags <- lapply(Xz,
                        function(x) Matrix::Diagonal(nrow(x), 1 / nrow(x)))
        I <- Matrix::bdiag(subdiags)
    } else {
        # all diagonal entries are 1
        n <- lapply(Xz, nrow) %>% Reduce(`+`, .)
        I <- Matrix::Diagonal(n)
    }
    return(I)
}

#' Create the q vector for an QP that solves min_x 0.5 * x'Px + q'x
#' @param kern_mats
#' @param target Vector of population means to re-weight to
#' @param aux_dim Dimension of auxiliary weights
#'
#' @return q vector
create_q_vector_multi_kern <- function(kern_mats, trtz) {
    J <- length(kern_mats)
    n <- Reduce(`+`, lapply(kern_mats, nrow))
    # concenate treated averages for each group
    q <- - do.call(c,
                 lapply(1:J,
                    function(j) colSums(kern_mats[[j]][trtz[[j]] == 1, , drop = F])
                       )
                )
    # q <- Matrix::sparseVector(q, 1:n,
    #                           2 * n)
    return(q)
}


#' Create the P matrix for an QP that solves min_x 0.5 * x'Px + q'x
#' @param X n x d matrix of covariates
#' @param Z Vector of group indicators
#'
#' @return P matrix
create_P_matrix_multi_kern <- function(kern_mats) {

    # up_right_mat <- Matrix::Matrix(c(0,1,0,0), ncol=2, byrow = T)

    # return(Matrix::kronecker(up_right_mat, Matrix::Diagonal(n)))
    return(Matrix::bdiag(kern_mats))
}


#' Create the constraints for QP: l <= Ax <= u
#' @param Xz list of J n x d matrices of covariates split by group
#' @param target Vector of population means to re-weight to
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#'
#' @return A, l, and u
create_constraints_multi_kern <- function(Xz, trtz, lowlim, uplim, 
                                          exact_global, verbose) {

    J <- length(Xz)
    # nj <- sapply(1:J, function(j) nrow(Xz[[j]]))
    n0j <- sapply(1:J, function(j) nrow(Xz[[j]]))
    n1j <- sapply(trtz, sum)
    d <- ncol(Xz[[1]])
    n <- Reduce(`+`, lapply(Xz, nrow))
    Xzt <- lapply(Xz, t)

    # dimension of auxiliary weights
    aux_dim <- n
  if(verbose) message("\tx Sum to one constraint")
  # sum-to-n1j constraint for each group
  A1 <- Matrix::t(Matrix::bdiag(lapply(Xz, function(x) rep(1, nrow(x)))))
  # A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = aux_dim))
  l1 <- n1j
  u1 <- n1j
  if(verbose) message("\tx Upper and lower bounds")
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  # A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
  l2 <- Reduce(c, sapply(1:J, function(j) rep(lowlim * n1j[[j]], n0j[j])))
  u2 <- Reduce(c, sapply(1:J, function(j) rep(uplim * n1j[[j]], n0j[j])))

  if(exact_global) {
      if(verbose) message("\tx Enforce exact global balance")
      # Constrain the overall mean to be equal to the target
      A3 <- do.call(cbind, lapply(1:J, function(j) Xzt[[j]]))
      # A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
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
    # if(verbose) message("\tx Sum to one constraint")
    # # sum-to-one constraint for each group
    # A1 <- Matrix::t(Matrix::bdiag(lapply(Xz, function(x) rep(1, nrow(x)))))
    # # A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = aux_dim))
    # l1 <- rep(1, J)
    # u1 <- rep(1, J)

    # if(verbose) message("\tx Upper and lower bounds")
    # # upper and lower bounds
    # A2 <- Matrix::Diagonal(n)
    # # A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = aux_dim))
    # l2 <- rep(lowlim, n)
    # u2 <- rep(uplim, n)

    # # exact_global <- FALSE
    # if(exact_global) {
    #     if(verbose) message("\tx Mantain overall population mean")
    #     # Constrain the overall mean to be equal to the target
    #     A3 <- do.call(cbind, lapply(1:J, function(j) Xzt[[j]] * n1j[j]))
    #     # A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = aux_dim))
    #     trt_sum <- Reduce(`+`,
    #         lapply(1:J,
    #             function(j) colSums(Xz[[j]][trtz[[j]] == 1, , drop = F])))
    #     l3 <- trt_sum
    #     u3 <- trt_sum
    # } else {
    #     if(verbose) message("\t(SKIPPING) Mantain overall population mean")
    #     # skip this constraint and just make empty
    #     A3 <- matrix(, nrow = 0, ncol = ncol(A2))
    #     l3 <- numeric(0)
    #     u3 <- numeric(0)
    # }

    # if(verbose) message("\tx Fit weights to data")
    # # constrain the auxiliary weights to be sqrt(P)'gamma
    # sqrtP <- Matrix::bdiag(Xzt)
    # A4 <- Matrix::cbind2(sqrtP, -Matrix::Diagonal(aux_dim))
    # l4 <- rep(0, aux_dim)
    # u4 <- rep(0, aux_dim)

    if(verbose) message("\tx Constrain treated weights to be zero")
    # zero out treated units
    A5 <- Matrix::bdiag(lapply(trtz, Matrix::diag))
    # A5 <- Matrix::cbind2(A5, Matrix::Matrix(0, nrow = nrow(A5), ncol = aux_dim))
    l5 <- numeric(n)
    u5 <- numeric(n)

    if(verbose) message("\tx Combining constraints")
    A <- rbind(A1, A2, A3, A5)
    l <- c(l1, l2, l3, l5)
    u <- c(u1, u2, u3, u5)

    return(list(A = A, l = l, u = u))
}
