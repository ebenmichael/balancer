

balancer <- function(X, trt, Z=NULL, V=NULL,
                     type=c("att", "subgrp", "subgrp_multi"),
                     link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                     regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc",
                                   "l1_all", "l1_nuc"),
                     hyperparam=NULL, nlambda=20, lambda.min.ratio=1e-3,
                     interact=F, normalized=TRUE,
                     ipw_weights=NULL, opts=list()) {
    #' Find Balancing weights by solving the dual optimization problem
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators or observed indicators
    #' @param V Group level covariates
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             subgroup ATTs with multilevel p-score, multilevel observational studies,
    #'             ATT with missing outcomes, and heterogeneous effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param hyperparam Regularization hyperparameter
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param interact Whether to interact group and individual level covariates
    #' @param Q m x m matrix to tie together ridge penalty, default: NULL,
    #'          if TRUE, use covariance of treated groups
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
        out <- balancer_att(X, trt, weightfunc, weightptr,
                            proxfunc, balancefunc, hyperparam,
                            nlambda, lambda.min.ratio,
                            ipw_weights, opts)
    } else if(type == "subgrp") {
        out <- balancer_subgrp(X, trt, Z, weightfunc, weightptr,
                                proxfunc, hyperparam, ridge, Q, NULL, NULL, opts)
    } else if(type == "subgrp_multi") {
        out <- balancer_multi(X, V, trt, Z, weightfunc, weightptr,
                              proxfunc, hyperparam, ridge, interact, opts)
    } else {
        stop("type must be one of ('att', 'subgrp', 'subgrp_multi')")
    }

    return(out)
}


balancer_att <- function(X, trt, weightfunc, weightfunc_ptr,
                         proxfunc, balancefunc, hyperparam=NULL,
                         nlambda=20, lambda.min.ratio=1e-3,
                         ipw_weights=NULL, opts=list()) {
    #' Balancing weights for ATT (in subgroups)
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param balancefunc Balance criterion measure
    #' @param hyperparam Regularization hyper parameter
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
                     ipw_weights=ipw_weights
                     )
    

    ## initialize at 0
    init = matrix(0, nrow=d, ncol=1)
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))
    
    
    ## if hyperparam is NULL, start from reference weights and decrease
    if(is.null(hyperparam)) {
        lam0 <- balancefunc(balancing_grad_att(init, loss_opts))
        lam1 <- lam0 * lambda.min.ratio
        ## decrease on log scale
        hyperparam <- exp(seq(log(lam0), log(lam1), length.out=nlambda))
    }


    ## collect results
    out <- list()
    out$theta <- matrix(,nrow=d, ncol=length(hyperparam))
    out$imbalance <- matrix(,nrow=d, ncol=length(hyperparam))    
    out$weights <- matrix(0, nrow=n, ncol=length(hyperparam))
    out$weightfunc <- weightfunc


    ## with multiple hyperparameters do warm starts        

    ## if(length(hyperparam) > 1) {

    prox_opts = list(lam=1)
    
    apgout <- apg_warmstart(make_balancing_grad_att(),
                            proxfunc, loss_opts, prox_opts,
                            hyperparam,
                            opts$x, opts$max_it, opts$eps,
                            opts$alpha, opts$beta, opts$accel, opts$verbose)

    ## weights and theta
    out$theta <- do.call(cbind, apgout)
    out$weights[trt==0,] <- apply(out$theta, 2,
                                  function(th) weightfunc(X[trt==0,,drop=F], as.matrix(th), ipw_weights))
    
    out$imbalance <- apply(out$theta, 2,
                           function(th) balancing_grad_att(as.matrix(th), loss_opts))

    out$lams <- hyperparam

    return(out)

}


balancer_subgrp <- function(X, trt, Z=NULL, weightfunc, weightfunc_ptr,
                            proxfunc, hyperparam, ridge, Q=NULL, kernel, kern_param=1,
                            ipw_weights=NULL, opts=list()) {
    #' Balancing weights for ATT (in subgroups)
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators    
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function, or "ridge"
    #' @param hyperparam Regularization hyper parameter
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

    ## if no subgroups, put everything into one group
    if(is.null(Z)) {
        Z <- rep(1, length(trt))
    }
    ## get the distinct group labels
    grps <- sort(unique(Z))
    m <- length(grps)

    n <- dim(X)[1]

    ## if(normalized) {
    ##     ## add a bias term
    ##     X <- cbind(rep(1, n), X)
    ## }

    d <- dim(X)[2]

    ## ##if no kernel use linear features
    ## if(is.null(kernel)) {
    ## get the group treated moments
    x_t <- sapply(grps,
                  function(k) colMeans(X[(trt ==1) & (Z==k), , drop=FALSE]))
    x_t <- as.matrix(x_t)
    
    Xc <- X[trt==0,,drop=FALSE]

    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     weight_type="subgroup",
                     z=Z[trt==0],
                     ridge=ridge,
                     ipw_weights=ipw_weights
                     )

    ## if ridge, make the identity matrix
    if(ridge) {
        loss_opts$hyper <- hyperparam
        if(!is.null(Q)) { 
            loss_opts$hasQ <- TRUE
            if(Q) {
                loss_opts$Q <- t(x_t)
            } else {
                loss_opts$Q <- Q
            }
        } else {
            loss_opts$hasQ <- FALSE
        }
    }
    
    ## indices of subgroups
    ctrl_z = Z[trt==0]
    loss_opts$z_ind <- lapply(grps, function(k) which(ctrl_z==k)-1)

    if(m==1) {
        loss_opts$weight_type="base"
    }

    if(length(hyperparam) == 1) {
        prox_opts = list(lam=hyperparam)
    } else if(length(hyperparam) == 2) {
        prox_opts = list(lam1=list(lam=hyperparam[1]),
                         lam2=list(lam=hyperparam[2]))
    } else {
        stop("hyperparam must be length 1 or 2")
    }

    ## initialize at 0
    init = matrix(0, nrow=d, ncol=m)

    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))
    if(is.null(kernel)) {
        apgout <- apg(make_balancing_grad(), proxfunc, loss_opts, prox_opts,
                      opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel, opts$verbose)
    } else {
        apgout <- apg(make_balancing_grad_kern(), proxfunc, loss_opts, prox_opts,
                      opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel, opts$verbose)
    }
                  

    ## collect results
    out <- list()

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    weights <- numeric(n)
    if(m == 1) {
        weights[trt==0] <- weightfunc(X[trt==0,,drop=F], theta, ipw_weights)
    } else {
        for(i in 1:m) {
            k = grps[i]
            weights[(trt == 0) & (Z == k)] <-
                weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE], theta[,i,drop=FALSE])
        }
    }
    out$weights <- weights

    ## The final imbalance
    loss_opts$ridge <- F
    out$imbalance <- balancing_grad(theta, loss_opts)
    out$weightfunc <- weightfunc
    return(out)
}





map_to_param <- function(link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                         regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc",
                                       "l1_all", "l1_nuc"),
                         ipw_weights=NULL,
                         normalized=F) {
    #' Map string choices to the proper parameters for the balancer sub-functions
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             subgroup ATTs with multilevel p-score, multilevel observational studies,
    #'             ATT with missing outcomes, and heterogeneous effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param normalized Whether to normalize the weights
    #'
    #' @return Parameters for balancer

    if(is.null(ipw_weights)) {
        ipw_weights = matrix(1, sum(trt==0), 1)    
    } else {
        ipw_weights = matrix(ipw_weights, sum(trt==0), 1)    
    }
    if(link == "logit") {
            if(normalized) {
                weightfunc <- exp_weights_ipw
                weightptr <- make_exp_weights_ipw()
            } else {
                weightfunc <- softmax_weights_ipw
                weightptr <- make_softmax_weights_ipw()
            }
    } else if(link == "linear") {
        weightfunc <- lin_weights_ipw
        weightptr <- make_lin_weights_ipw()
    } else if(link == "pos-linear") {
        weightfunc <- pos_lin_weights_ipw
        weightptr <- make_pos_lin_weights_ipw()        
    } else if(link == "enet") {
        stop("Elastic Net not impleneted")
    } else if(link == "pos-enet") {
        stop("Elastic Net not impleneted")
    } else {
        stop("link must be one of ('logit', 'linear', 'pos-linear')")
    }

    if(is.null(regularizer)) {
        proxfunc <- make_no_prox()
    } else if(regularizer == "l1") {
        proxfunc <- make_prox_l1()
        balancefunc <- linf
    } else if(regularizer == "grpl1") {
        proxfunc <- make_prox_l1_grp()
    } else if(regularizer == "l1grpl1") {
        proxfunc <- make_prox_l1_grp_l1()
        ## double the covariate matrix to include two sets of parameters
        X <- cbind(X,X)
    } else if(regularizer == "l2") {
        proxfunc <- make_prox_l2()
        balancefunc <- l2
    } else if(regularizer == "ridge") {
        proxfunc <- make_prox_l2_sq()
        balancefunc <- l2sq
    } else if(regularizer == "linf") {
        stop("L infinity regularization not implemented")
    } else if(regularizer == "nuc") {
        proxfunc <- make_prox_nuc()
    } else if(regularizer == "nucl1") {
        proxfunc <- make_prox_nuc_l1()
        ## double the covariate matrix to include two sets of parameters
        X <- cbind(X,X)
    } else if(regularizer == "l1_all") {
        proxfunc <- make_prox_l1_all()
    } else if(regularizer == "l1_nuc") {
        proxfunc <- make_prox_l1_nuc()
    } else {
        stop("regularizer must be one of (NULL, 'l1', 'grpl1', 'l1grpl1', 'ridge', 'linf', 'nuc', 'nucl1')")
    }


    return(list(weightfunc, weightptr, proxfunc, balancefunc, ipw_weights))
}
