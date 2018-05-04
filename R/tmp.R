
logsumexp <- function(x0) {
    #' Compute numerically stable logsumexp
    m <- max(x0)
    val <- log(sum(exp(x0 - m))) + m
    return(val)
}


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

balancer_subgrp2 <- function(X, trt, Z=NULL, weightfunc,
                            proxfunc, hyperparam, opts=list()) {
    #' Helper function to fit the dual for general odds function and prox
    #' estimates heterogeneous treatment effects
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
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
    
    ## objective gradient
    ## f^*'(x * theta) - x
    grad <- function(theta) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)
        ## first part of gradient comes from control units
        grad1 <- sapply(grps,
                        function(k) t(X[(trt == 0) & (Z == k),]) %*%
                                    weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE] %*% theta[,k]))

        ## second part of gradient comes from treated units
        grad <- grad1 - x_t

        return(as.numeric(grad))
    }

    ## prox operator of regularizer

    prox <- function(theta, step) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)

        ## apply prox operator for covariate parameters
        proxtheta <- proxfunc(theta, step * hyperparam)
        
        return(as.numeric(proxtheta))
    }

    loss <- function(theta) {
        l <- logsumexp(X[trt==0,] %*% theta) - t(x_t) %*% theta
        ## l <- sum((X[trt==0,] %*% theta)^2) - t(x_t) %*% theta
        ## print(as.numeric(l))
        return(as.numeric(l))
    }

    ## combine opts with defaults
    ## opts <- c(max_it=5000, eps=1e-8, beta=.9, opts)
    print(opts$eps)
    print(opts$max_it)

    ## apgout <- apg2(loss, grad, prox, d * m, opts$max_it, opts$eps, opts$beta, opts$accel)
    apgout <- apg3(grad, prox, d * m, opts$max_it, opts$eps, opts$beta, opts$accel, opts$alpha)
    ## collect results
    out <- list()

    ## theta
    theta <- matrix(apgout, ncol=m)
    out$theta <- theta
    ## weights
    weights <- numeric(n)
    for(k in grps) {
        weights[(trt == 0) & (Z == k)] <-
            weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE] %*% theta[,k])
    }

    out$weights <- weights

    ## The final imbalance
    out$imbalance <- matrix(grad(theta), ncol=m)
    return(out)
}
