

balancer2 <- function(X, trt, Z=NULL, type=c("att", "subgrp", "missing", "hte"),
                     link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                     regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc"),
                     hyperparam, normalized=TRUE, opts=list()) {
    #' Find Balancing weights by solving the dual optimization problem
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators or observed indicators
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             ATT with missing outcomes, and heterogeneouts effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param hyperparam Regularization hyperparameter
    #' @param normalized Whether to normalize the weights
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

    if(link == "logit") {
        weightfunc <- softmax
    } else if(link == "linear") {
        if(normalized) {
            weightfunc <- normlin
        } else {
            weightfunc <- identity
        }
    } else if(link == "pos-linear") {
        if(normalized) {
            weightfunc <- normposlin
        } else {
            weightfunc <- poslin
        }
    } else if(link == "enet") {
        alpha <- if(is.null(opts$alpha)) 0.5 else opts$alpha
        weightfunc <- function(eta) enet(eta, alpha)
    } else if(link == "pos-enet") {
        alpha <- if(is.null(opts$alpha)) 0.5 else opts$alpha
        weightfunc <- function(eta) posenet(eta, alpha)
    } else {
        stop("link must be one of ('logit', 'linear', 'pos-linear')")
    }

    if(is.null(regularizer)) {
        proxfunc <- no_prox
    } else if(regularizer == "l1") {
        proxfunc <- l1_prox
    } else if(regularizer == "grpl1") {
        proxfunc <- l1_grp_prox
    } else if(regularizer == "l2") {
        proxfunc <- l2_prox
    } else if(regularizer == "ridge") {
        proxfunc <- ridge_prox
    } else if(regularizer == "linf") {
        proxfunc <- linf_prox
    } else if(regularizer == "nuc") {
        proxfunc <- nuc_prox
    } else {
        stop("regularizer must be one of (NULL, ;l1', 'grpl1', 'ridge', 'linf', 'nuc')")
    }
    if(type == "att") {
        out <- balancer_subgrp2(X, trt, NULL, weightfunc, proxfunc, hyperparam, opts)
    } else if(type == "subgrp") {
        out <- balancer_subgrp2(X, trt, Z, weightfunc, proxfunc, hyperparam, opts)
    } else if(type == "missing") {
        out <- balancer_missing(X, trt, Z, weightfunc, proxfunc, hyperparam, opts)
    } else if(type == "hte") {
        out <- balancer_hte(X, trt, weightfunc, proxfunc, hyperparam, opts)
    } else {
        stop("type must be one of ('att', 'subgrp', 'missing', 'hte')")
    }

    return(out)
}



balancer_subgrp2 <- function(X, trt, Z=NULL, weightfunc, weightfunc_ptr,
                             proxfunc, hyperparam, opts=list()) {
    #' Balancing weights for ATT (in subgroups)
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators    
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param hyperparam Regularization hyper parameter
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

    ## get the group treated moments
    x_t <- sapply(grps,
                    function(k) colMeans(X[(trt ==1) & (Z==k), , drop=FALSE]))

    ##x_t <- as.numeric(x_t)
    
    loss_opts = list(Xc=X[trt==0,],
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     weight_type="subgroup",
                     z=Z[trt==0]
                     )
    if(m==1) {
        loss_opts$weight_type="base"
    }

    prox_opts = list(lam=hyperparam)

    ## initialize at 0
    init = matrix(0, nrow=dim(X)[2], ncol=m)
    
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=50000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init))
    

    apgout <- apg(make_balancing_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel)
                  

    ## collect results
    out <- list()

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    weights <- numeric(n)
    for(i in 1:m) {
        k = grps[i]
        weights[(trt == 0) & (Z == k)] <-
            weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE], theta[,i,drop=FALSE])
    }

    out$weights <- weights

    ## The final imbalance
    out$imbalance <- balancing_grad(theta, loss_opts)
    return(out)
}



balancer_missing2 <- function(X, trt, R, weightfunc, weightfunc_ptr,
                             proxfunc, hyperparam, opts=list()) {
    #' Balancing weights for missing outcomes
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param R vector of missingness indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param hyperparam Regularization hyper parameter
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}}
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}

    m <- 2

    n <- dim(X)[1]

    d <- dim(X)[2]

    ## get the moments for treated units twice
    x_t <- cbind(colMeans(X[trt==1,]),
                 colMeans(X[trt==1,]))
    loss_opts = list(Xc=X[R==1,],
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     weight_type="missing",
                     trt=trt[R==1]
                     )
    prox_opts = list(lam=hyperparam)


    ## initialize at 0
    init = matrix(0, nrow=dim(X)[2], ncol=m)
    
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=50000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init))
    

    apgout <- apg(make_balancing_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel)
                  

    ## collect results
    out <- list()

    ## theta
    theta <- apgout
    out$theta <- theta

    ## weights    
    ## weights for R=1 T=0 to T=1
    weights1 <- numeric(n)
    weights1[trt == 0 & R == 1] <-
        weightfunc(X[trt == 0 & R == 1,,drop=FALSE], theta[,1, drop=FALSE])
    
    ## weights for R=1 T=1 to T=1
    weights2 <- numeric(n)
    weights2[trt == 1 & R == 1] <-
        weightfunc(X[trt == 1 & R == 1,,drop=FALSE], theta[,2, drop=FALSE])
    weights <- cbind(weights1, weights2)

    out$weights <- weights

    ## The final imbalance
    out$imbalance <- balancing_grad(theta, loss_opts)
    return(out)
}




balancer_hte2 <- function(X, trt, weightfunc, weightfunc_ptr,
                             proxfunc, hyperparam, opts=list()) {
    #' Balancing weights for heterogeneous treatment effects
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param hyperparam Regularization hyper parameter
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}}
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}
    
    ## one set of weights for each treated unit
    m <- sum(trt)

    n <- dim(X)[1]
    
    d <- dim(X)[2]
    
    ## keep the covariates for the treated units
    x_t <- t(X[trt==1,])

    ## initialize at 0
    init = matrix(0, nrow=dim(X)[2], ncol=m)
    
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=50000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init))
    
    loss_opts = list(Xc=X[trt==0,],
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     weight_type="hte"
                     )
    prox_opts = list(lam=hyperparam)

    apgout <- apg(make_balancing_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel)
                  

    ## collect results
    out <- list()

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    eta <- X %*% theta
    eta[trt==1,] <- 0
    weights <- weightfunc(eta)
    out$weights <- weights

    ## The final imbalance
    out$imbalance <- balancing_grad(theta, loss_opts)
    return(out)
}
