######################################################################
## Approximate calibration a la SBW
######################################################################

#' Approximate calibration with categorical covariates
#' 
calibrate_approx <- function(X, target, lambda, lowlim = 0, uplim = 1,
                             verbose = TRUE, eps_abs = 1e-5, eps_rel = 1e-5, 
                             ...) {
    n <- nrow(X)
    d <- ncol(X)
    # Setup the components of the QP and solve
    # q vector is zero
    q_vec <- Matrix::sparseVector(0, 0, n)
    # P matrix is identity
    P_mat <- Matrix::.symDiagonal(n)

    if(verbose) message("Creating constraint matrix...")
    constraints <- create_calibrate_constraints(X, target, lambda, lowlim,
                                                uplim, verbose)
    
    # solve the QP
    settings <- do.call(osqp::osqpSettings, 
                        c(list(verbose = verbose, 
                               eps_rel = eps_rel,
                               eps_abs = eps_abs), 
                        list(...)))
    solution <- osqp::solve_osqp(P_mat, q_vec, constraints$A,
                                constraints$l, constraints$u,
                                pars = settings)
    # primal weights and dual parameters
    weights <- solution$x
    dual <- c(solution$y[1],
              solution$y[(length(solution$y) - ncol(X) + 1):length(solution$y)])
    # balance
    imbalance <- target - Matrix::t(X) %*% weights
    return(list(weights = weights, 
                dual = dual,
                imbalance = imbalance))



}


#' Create the constraints for QP: l <= Ax <= u
#' @param X Matrix of factors'
#' @param target Vector of population means to re-weight to
#' @param lambda Balance constraint level
#' @param lowlim Lower limit on weights
#' @param uplim Upper limit on weights
#'
#' @return A, l, and u
create_calibrate_constraints <- function(X, target, lambda, 
                                        lowlim, uplim, verbose) {

    n <- nrow(X)
    n_fac <- Matrix::colSums(X)

    if(verbose) message("\tx Sum to one constraint")
    # sum-to-one constraint for each group
    A1 <- Matrix::Matrix(1, nrow = 1, ncol = n)
    l1 <- 1
    u1 <- 1

    if(verbose) message("\tx Upper and lower bounds")
    # upper and lower bounds
    A2 <- Matrix::Diagonal(n)
    l2 <- rep(lowlim, n)
    u2 <- rep(uplim, n)

    if(verbose) message("\tx Balance constraints")
    A3 <- Matrix::t(X)
    # element wise max with 0
    l3 <- pmax(target - lambda / n_fac, 0)
    u3 <- target + lambda / n_fac


    if(verbose) message("\tx Combining constraints")
    A <- rbind(A1, A2, A3)
    l <- c(l1, l2, l3)
    u <- c(u1, u2, u3)

    return(list(A = A, l = l, u = u))
}
