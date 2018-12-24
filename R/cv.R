################################################################################
## Cross validation
################################################################################

create_folds <- function(n, k) {
    #' Partition indices into k folds
    #' @param n Number of data points
    #' @param k Number of folds
    #'
    #' @return List of fold indices

    ## permute
    perm <- sample(n)
    ## partition
    return(lapply(1:k, function(i) perm[seq(i, n, by=k)]))
}


balancer_cv <- function(X, trt, y=NULL, k=10, Z=NULL, V=NULL,
                        type=c("att", "subgrp", "subgrp_multi"),
                        link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                        regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc",
                                      "l1_all", "l1_nuc"),
                        lambda=NULL, nlambda=20, lambda.min.ratio=1e-3,
                        interact=F, normalized=TRUE,
                        ipw_weights=NULL, opts=list()) {
    #' Find Balancing weights by solving the dual optimization problem
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param y Vector of outcomes to estimate effect(s). If NULL then only return weights
    #' @param k Number of folds
    #' @param Z Vector of subgroup indicators or observed indicators
    #' @param V Group level covariates
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             subgroup ATTs with multilevel p-score, multilevel observational studies,
    #'             ATT with missing outcomes, and heterogeneous effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param lambda Vector of regularization hyperparameters to consider, if NULL then uses
    #'               a data-generated list
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param interact Whether to interact group and individual level covariates
    #' @param normalized Whether to normalize the weights
    #' @param ipw_weights Separately estimated IPW weights to measure dispersion against, default is NULL
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
    params <- map_to_param(link, regularizer, ipw_weights, normalized)
    weightfunc <- params[[1]]
    weightptr <- params[[2]]
    proxfunc <- params[[3]]
    balancefunc <- params[[4]]
    ipw_weights <- params[[5]]
    if(type == "att") {
        out <- balancer_att_cv(X, trt, k, y, weightfunc, weightptr,
                            proxfunc, balancefunc, lambda,
                            nlambda, lambda.min.ratio,
                            ipw_weights, opts)
    } else {
        stop("type must be one of ('att')")
    }

    return(out)
}


balancer_att_cv <- function(X, trt, k, y=NULL, weightfunc, weightfunc_ptr,
                            proxfunc, balancefunc, lambda=NULL,
                            nlambda=20, lambda.min.ratio=1e-3,
                            ipw_weights=NULL, opts=list()) {
    #' Balancing weights for ATT (in subgroups)
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param k Number of folds
    #' @param y Vector of outcomes to estimate effect(s). If NULL then only return weights    
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param balancefunc Balance criterion measure
    #' @param lambda Vector of regularization hyperparameters to consider, if NULL then uses
    #'               a data-generated list
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param ridge Whether to use L2 penalty in dual
    #' @param Q m x m matrix to tie together ridge penalty, default: NULL,
    #'          if TRUE, use covariance of treated groups
    #' @param kernel What kernel to use, default NULL
    #' @param kern_param Hyperparameter for kernel
    #' @param ipw_weights Separately estimated IPW weights to measure dispersion against, default is NULL
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}}
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}

    n <- dim(X)[1]

    d <- dim(X)[2]

    x_t <- matrix(colSums(X[trt==1,,drop=FALSE]), nrow=d)
    Xc <- X[trt==0,,drop=FALSE]
    
    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     ipw_weights=ipw_weights[trt==0,,drop=F]
                     )

    ## initialize at 0
    init = matrix(0, nrow=d, ncol=1)    
    
    ## if hyperparam is NULL, start from reference weights and decrease
    if(is.null(lambda)) {
        lam0 <- balancefunc(balancing_grad_att(init, loss_opts))
        lam1 <- lam0 * lambda.min.ratio
        ## decrease on log scale
        lambda <- exp(seq(log(lam0), log(lam1), length.out=nlambda))
    }


    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))

    prox_opts = list(lam=1)

    folds <- create_folds(n, k)

    ## helper function to fit weights and get balance
    fit_fold <- function(i) {

        ## subset data 
        idx <- folds[[i]]
        
        loss_opts$Xc <- X[-idx,][trt[-idx]==0,,drop=FALSE]
        loss_opts$Xt <- matrix(colSums(X[-idx,][trt[-idx]==1,,drop=FALSE]), nrow=d)
        loss_opts$ipw_weights <- ipw_weights[-idx,,drop=F][trt[-idx]==0,,drop=F]

        ## fit balancing weights
        apgout <- apg_warmstart(make_balancing_grad_att(),
                                proxfunc, loss_opts, prox_opts,
                                lambda,
                                opts$x, opts$max_it, opts$eps,
                                opts$alpha, opts$beta, opts$accel, opts$verbose)

        ## evaluate imbalance on held out data
        loss_opts$Xc <- X[idx,][trt[idx]==0,,drop=FALSE]
        loss_opts$Xt <- matrix(colSums(X[idx,][trt[idx]==1,,drop=FALSE]), nrow=d)
        loss_opts$ipw_weights <- ipw_weights[idx,,drop=F][trt[idx]==0,,drop=F]

        vapply(1:length(lambda),
               function(i) balancefunc(balancing_grad_att(apgout[[i]], loss_opts)),
               numeric(1))
    }

    ## fit and evaluate
    bals <- vapply(1:k, fit_fold, numeric(length(lambda)))

    ## fit on whole data
    ipw_weights <- ipw_weights[trt==0,,drop=F]
    loss_opts$Xc <- X[trt==0,,drop=FALSE]
    loss_opts$Xt <- matrix(colSums(X[trt==1,,drop=FALSE]), nrow=d)
    loss_opts$ipw_weights <- ipw_weights

    ## fit balancing weights
    apgout <- apg_warmstart(make_balancing_grad_att(),
                            proxfunc, loss_opts, prox_opts,
                            lambda,
                            opts$x, opts$max_it, opts$eps,
                            opts$alpha, opts$beta, opts$accel, opts$verbose)
    
    
    ## collect results
    out <- list()
    out$theta <- matrix(,nrow=d, ncol=length(lambda))
    out$imbalance <- matrix(,nrow=d, ncol=length(lambda))    
    out$weights <- matrix(0, nrow=n, ncol=length(lambda))
    out$weightfunc <- weightfunc



    ## weights and theta
    out$theta <- do.call(cbind, apgout)
    out$weights[trt==0,] <- apply(out$theta, 2,
                                  function(th) weightfunc(X[trt==0,,drop=F], as.matrix(th), ipw_weights))
    
    out$imbalance <- apply(out$theta, 2,
                           function(th) balancing_grad_att(as.matrix(th), loss_opts))

    out$lambda <- lambda

    out$cvm <- rowMeans(bals)
    out$cvsd <- apply(bals, 1, sd)
    out$cv <- bals

    ## estimate treatment effect(s)
    out$att <- mean(y[trt==1]) - t(y) %*% out$weights / sum(trt==1)

    
    return(out)

}
