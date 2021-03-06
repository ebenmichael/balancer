################################################################################
## Calibrate multiway tables
################################################################################

#' Calibrate sample to target
#' @param data Dataframe with covariate information, sample and target counts
#' @param formula sample_count ~ covariates
#' @param target_count Number in cell in target
#' @param order What order interactions to balance
#' @param lambda Regularization hyperparamter
#' @param lowlim Lower bound on weights, default 0
#' @param uplim Upper bound on weights, default Inf
#' @param prob_weights Optional sampling weights to include
#' @export
calibrate <- function(formula, target_count, data,
                      order = NULL, lambda = 1,
                      lambda_max = NULL, n_lambda = 100,
                      lowlim = 0, uplim = Inf, 
                      prob_weights = NULL,
                      verbose = FALSE, ...) {
  
  # create distinct cells for all interactions
  if(verbose) message("Creating table of cell counts")
  cells <- create_cells(formula, enquo(target_count), data)

  # get weights
  weights <- calibrate_(cells %>% select(-sample_count, -target_count),
                        cells$sample_count, cells$target_count,
                        order, lambda, lambda_max, n_lambda,
                        lowlim, uplim, prob_weights, verbose, ...)
  # combine back in and return
  cells %>% filter(sample_count != 0) %>%
    bind_cols(weights) %>%
    select(-sample_count, -target_count) %>%
    right_join(cells,
               by = cells %>% select(-sample_count, -target_count) %>%
                    names()) %>%
    pivot_longer(!names(cells), names_to = "lambda", values_to = "weight") %>%
    mutate(weight = replace_na(weight, 0),
           lambda = as.numeric(lambda)) %>%
    return()
}

create_cells <- function(formula, target_count, data) {


  covs <- all.vars(terms(Formula::Formula(formula), rhs = 1)[[3]])
  sample_count <- terms(Formula::Formula(formula), rhs = 1)[[2]]
  cells <- data %>%
    mutate(across(covs, as.factor)) %>%
    group_by(across(covs)) %>%
    summarise(sample_count = sum(!!sample_count),
              target_count = sum(!!target_count)) %>%
    ungroup()

  return(cells)
  # cells <- cells %>% select(-sample_count, -target_count)
  # return(list(cells = cells, sample_counts = sample_counts,
  #             target_counts = target_counts))
}

calibrate_ <- function(cells, sample_counts, target_counts, order = NULL,
                      lambda = 1, lambda_max = NULL, n_lambda = 100,
                      lowlim = 0, uplim = Inf,
                      prob_weights = NULL, verbose = FALSE,
                      ...) {

  if(verbose) message("Creating design matrix")
  # get design matrix for the number of interactions
  D <- create_design_matrix(cells, order)

  # create constraints for raking
  constraints <- create_rake_constraints(cells, D, sample_counts,
                                         target_counts, lowlim, uplim,
                                         verbose)

  # return(constraints)
  # P matrix and q vector
  if(verbose) message("Creating quadratic term matrix")
  if(is.null(lambda)) {
    if(is.null(lambda_max)) {
      unif_imbal <- Matrix::t(D) %*% (sample_counts - target_counts)
      lambda_max <- sqrt(sum(unif_imbal ^ 2))
    }
    lam_seq <- lambda_max * 10 ^ seq(0, -5, length.out = n_lambda)
    P <- create_rake_Pmat(D, sample_counts, 0)
  } else {
    P <- create_rake_Pmat(D, sample_counts, lambda)
  }
  
  
  if(verbose) message("Creating linear term vector")
  qvec <- create_rake_qvec(D, sample_counts, target_counts, lambda,
                           prob_weights)


  settings <- do.call(osqp::osqpSettings,
                        c(list(...), list(verbose = verbose)))
  if(verbose) message("Setting up optimization")
  solver <- osqp::osqp(P, qvec, constraints$A, constraints$l, constraints$u,
                       pars = settings)
  if(verbose) message("Optimizing")
  # solution <- osqp::solve_osqp(P, qvec, constraints$A,
  #                              constraints$l, constraints$u,
  #                              pars = settings)
  if(is.null(lambda)) {
    
    # get diagonal elements of P
    # P_diag <- compute_Pmat_diag(D, sample_counts)
    # diag_idx <- sapply(1: length(P_diag), function(i) c(i,i))
    # print(diag_idx[1:10])
      
    x <- NULL
    y <- NULL
    solution <- matrix(, nrow = nrow(P), ncol = n_lambda)
    for(i in 1:length(lam_seq)) {
      if(verbose) message(paste0("Solving with lambda = ", lam_seq[i]))
      
      P_new <- update_rake_Pmat(P, sample_counts, lam_seq[i])
      solver$Update(Px = Matrix::triu(P_new)@x)
      # solver$Update(Px = Matrix::diag(P_new), Px_idx = 1:(nrow(P_new)))
      solver$WarmStart(x = x, y = y)
      sol_lam <- solver$Solve()
      x <- sol_lam$x
      y <- sol_lam$y
      solution[, i] <- x
    }
    solution <- as.data.frame(solution)
    names(solution) <- lam_seq
  } else {
    solution <- data.frame(solver$Solve()$x)
    names(solution) <- lambda
  }
  return(solution)

}

create_design_matrix <- function(cells, order) {

  if(is.null(order)) {
    order <- ncol(cells)
  }
  if(order == 1) {
    return(Matrix::sparseMatrix(dims = c(nrow(cells), 1), i = {}, j = {}))
  }

  contrast_cells <- lapply(cells, contrasts, contrasts = F)

  form <- as.formula(paste("~ . ^ ", order, " - . + 0"))

  D <- Matrix::sparse.model.matrix(form, data = cells)
  return(D)
}


create_arrays <- function(cells, sample_counts, target_counts) {

  dnames <- lapply(cells, levels)
  dims <- sapply(dnames, length)

  sample_arr <- array(NA, dim = dims, dimnames = dnames)
  target_arr <- array(NA, dim = dims, dimnames = dnames)

  sample_arr[as.matrix(cells)] <- sample_counts
  target_arr[as.matrix(cells)] <- target_counts

  return(list(
        sample = tidyr::replace_na(sample_arr, 0),
        target = tidyr::replace_na(target_arr, 0)
        ))
}


create_rake_constraints <- function(cells, D, sample_counts, target_counts,
                                    lowlim, uplim, verbose) {

  if(verbose) message("Creating constraint matrix")
  # number of cells
  n_cells <- nrow(cells)
  # number of non-empty cells
  n_nempty <- sum(sample_counts != 0)
  sample_counts_nempty <- sample_counts[sample_counts != 0]

  # sum to N constraint
  if(verbose) message("\tx Sum to N constraint")
  # A_sumN <- Matrix::cbind2(
  #     Matrix::Matrix(sample_counts, nrow = 1, ncol = n_cells),
  #     Matrix::Matrix(0, nrow = 1, ncol = ncol(D))
  #     )
  A_sumN <- Matrix::Matrix(sample_counts_nempty, nrow = 1, ncol = n_nempty)
  l_sumN <- sum(target_counts)
  u_sumN <- l_sumN


  # non-negative constraint
  if(verbose) message("\tx Non-negativity constraint")
  # A_nn <- Matrix::cbind2(Matrix::Diagonal(n_cells),
  #                       Matrix::Matrix(0, nrow = n_cells, ncol = ncol(D)))
  A_nn <- Matrix::Diagonal(n_nempty)
  # l_nn <- rep(0, n_cells)
  # u_nn <- rep(Inf, n_cells)
  l_nn <- rep(lowlim, n_nempty)
  u_nn <- rep(uplim, n_nempty)

  # marginal constraints
  if(verbose) message("\tx Exact marginal balance constraints")
  design_mat <- Matrix::sparse.model.matrix(~ . - 1, cells)
  # A_marg <- Matrix::cbind2(Matrix::t(design_mat * sample_counts),
  #                         Matrix::Matrix(0, nrow = ncol(design_mat),
  #                                        ncol = ncol(D)))
  A_marg <- Matrix::t(design_mat[sample_counts != 0, , drop = F] * 
                      sample_counts_nempty)
  l_marg <- as.numeric(Matrix::t(design_mat) %*% target_counts)
  u_marg <- l_marg

  # pre-compute t(D) %*% w
  # if(verbose) message("\tx Pre-computing design matrix multiplication")
  # A_D <- Matrix::cbind2(Matrix::t(D * sample_counts), -Matrix::Diagonal(ncol(D)))
  # l_D <- numeric(ncol(D))
  # u_D <- numeric(ncol(D))

  # combine
  A <- rbind(A_sumN, A_nn, A_marg)#, A_D)
  l <- c(l_sumN, l_nn, l_marg)#, l_D)
  u <- c(u_sumN, u_nn, u_marg)#, u_D)

  return(list(A = A, l = l, u = u))

}

create_rake_Pmat <- function(D, sample_counts, lambda) {

  # P <- Matrix::bdiag(c(as.list(rep(lambda, nrow(D))),
  #                      as.list(rep(1, ncol(D)))))
  nnz_idxs <- which(sample_counts != 0)
  sqrtP <- D[nnz_idxs,, drop = F] * sample_counts[nnz_idxs]
  P <- sqrtP %*% Matrix::t(sqrtP) +
    lambda * sample_counts[nnz_idxs] * Matrix::Diagonal(length(nnz_idxs))
  return(P)

}

update_rake_Pmat <- function(P, sample_counts, lambda) {

  nnz_idxs <- which(sample_counts != 0)
  P <- P +
    lambda * sample_counts[nnz_idxs] * Matrix::Diagonal(length(nnz_idxs))
  return(P)
}



compute_Pmat_diag <- function(D, sample_counts) {

  nnz_idxs <- which(sample_counts != 0)
  d_counts <- Matrix::diag(D[nnz_idxs,, drop = F] %*% Matrix::t(D[nnz_idxs,, drop = F]))
  return(d_counts * sample_counts[nnz_idxs] ^ 2)
}

create_rake_qvec <- function(D, sample_counts, target_counts, lambda, prob_weights) {

  nnz_idxs <- which(sample_counts != 0)

  if(is.null(prob_weights)) {
    # prob_weights <- numeric(nrow(D))
    prob_weights <- numeric(length(nnz_idxs))
  }

  rhs <- Matrix::t(D) %*% target_counts
  # lhs <- sample_counts * D
  lhs <- D[nnz_idxs,, drop = F] * sample_counts[nnz_idxs]

  # return(-c(as.numeric(lhs %*% rhs) + lambda * prob_weights,
  #           numeric(ncol(D))))
  return(-c(as.numeric(lhs %*% rhs)))# + lambda * prob_weights))
}