################################################################################
## Calibrate multiway tables
################################################################################

#' Calibrate sample to target
#' @param data Dataframe with covariate information, sample and target counts
#' @param formula sample_count ~ covariates
#' @param target_count Number in cell in target
#' @param order What order interactions to balance
#' @param lambda Regularization hyperparamter
#' @param prob_weights Optional sampling weights to include
#' @export
calibrate <- function(formula, target_count, data,
                      order = NULL, lambda = 1, prob_weights = NULL,
                      verbose = FALSE, ...) {
  
  # create distinct cells for all interactions
  if(verbose) message("Creating table of cell counts")
  cells <- create_cells(formula, enquo(target_count), data)

  # get weights
  weights <- calibrate_(cells %>% select(-sample_count, -target_count),
                        cells$sample_count, cells$target_count,
                        order, lambda, prob_weights, verbose, ...)

  # combine back in and return
  cells %>% filter(sample_count != 0) %>%
    mutate(weight = weights) %>%
    select(-sample_count, -target_count) %>%
    right_join(cells,
               by = cells %>% select(-sample_count, -target_count) %>%
                    names()) %>%
    mutate(weight = replace_na(weight, 0)) %>%
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
                      lambda = 1, prob_weights = NULL, verbose = FALSE,
                      ...) {

  if(verbose) message("Creating design matrix")
  # get design matrix for the number of interactions
  D <- create_design_matrix(cells, order)

  # create constraints for raking
  constraints <- create_rake_constraints(cells, D, sample_counts, 
                                         target_counts, verbose)

  # return(constraints)
  # P matrix and q vector
  if(verbose) message("Creating quadratic term matrix")
  P <- create_rake_Pmat(D, sample_counts, lambda)
  
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
  solution <- solver$Solve()
  return(solution$x)
  # return(solution$x[1:nrow(cells)])

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
                                    verbose) {

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
  l_nn <- rep(0, n_nempty)
  u_nn <- rep(Inf, n_nempty)

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
  P <- sqrtP %*% Matrix::t(sqrtP) + lambda * Matrix::Diagonal(length(nnz_idxs))
  return(P)

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
  return(-c(as.numeric(lhs %*% rhs) + lambda * prob_weights))
}