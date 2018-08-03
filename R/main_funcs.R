

balancer <- function(X, trt, Z=NULL, V=NULL,
                     type=c("att", "subgrp", "subgrp_multi",
                            "multimatch", "missing", "hte"),
                     link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                     regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc",
                                   "l1_all", "l1_nuc"),
                     hyperparam, interact=F, Q=NULL, kernel=NULL, kern_param=1, normalized=TRUE, opts=list()) {
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
    #' @param interact Whether to interact group and individual level covariates
    #' @param Q m x m matrix to tie together ridge penalty, default: NULL,
    #'          if TRUE, use covariance of treated groups
    #' @param kernel What kernel to use, default NULL
    #' @param kern_param Hyperparameter for kernel
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



    if(type != "hte") {
        if(link == "logit") {
            if(type %in% c("subgrp_multi", "multi_match")) {
                weightfunc <- exp_weights
                weightptr <- make_exp_weights()
            }
            else {
                if(normalized) {
                    weightfunc <- softmax_weights
                    weightptr <- make_softmax_weights()
                } else {
                    weightfunc <- softmax_weights
                    weightptr <- make_softmax_weights()
                }
            }
        } else if(link == "linear") {
            if(normalized) {
                weightfunc <- lin_weights
                weightptr <- make_lin_weights()
            } else {
                weightfunc <- lin_weights
                weightptr <- make_lin_weights()
            }
        } else if(link == "pos-linear") {
            if(normalized) {
                weightfunc <- pos_lin_weights
                weightptr <- make_pos_lin_weights()
            } else {
                weightfunc <- pos_lin_weights
                weightptr <- make_pos_lin_weights()
            }
        } else if(link == "enet") {
            stop("Elastic Net not impleneted")
        } else if(link == "pos-enet") {
            stop("Elastic Net not impleneted")
        } else {
            stop("link must be one of ('logit', 'linear', 'pos-linear')")
        }
    } else {
        linear <- FALSE
        if(link == "logit") {
            if(type %in% c("subgrp_multi", "multi_match")) {
                weightfunc <- exp_weights2
                weightptr <- make_exp_weights2()
            }
            else {
                if(normalized) {
                    weightfunc <- softmax_weights2
                    weightptr <- make_softmax_weights2()
                } else {
                    weightfunc <- softmax_weights2
                    weightptr <- make_softmax_weights2()
                }
            }
        } else if(link == "linear") {
            if(normalized) {
                weightfunc <- lin_weights2
                weightptr <- make_lin_weights2()
            } else {
                weightfunc <- lin_weights2
                weightptr <- make_lin_weights2()
            }
            linear <- TRUE
        } else if(link == "pos-linear") {
            if(normalized) {
                weightfunc <- pos_lin_weights2
                weightptr <- make_pos_lin_weights2()
            } else {
                weightfunc <- pos_lin_weights2
                weightptr <- make_pos_lin_weights2()
            }
        } else if(link == "enet") {
            stop("Elastic Net not impleneted")
        } else if(link == "pos-enet") {
            stop("Elastic Net not impleneted")
        } else {
            stop("link must be one of ('logit', 'linear', 'pos-linear')")
        }        
    }

    ridge <- FALSE
    if(is.null(regularizer)) {
        proxfunc <- make_no_prox()
    } else if(regularizer == "l1") {
        proxfunc <- make_prox_l1()
    } else if(regularizer == "grpl1") {
        proxfunc <- make_prox_l1_grp()
    } else if(regularizer == "l1grpl1") {
        proxfunc <- make_prox_l1_grp_l1()
        ## double the covariate matrix to include two sets of parameters
        X <- cbind(X,X)
    } else if(regularizer == "l2") {
        proxfunc <- make_prox_l2()
    } else if(regularizer == "ridge") {
        proxfunc <- make_no_prox()
        ridge <- TRUE
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
    
    if(type == "att") {
        out <- balancer_subgrp(X, trt, NULL, weightfunc, weightptr,
                               proxfunc, hyperparam, ridge, Q, kernel,
                               kern_param, opts)
    } else if(type == "subgrp") {
        out <- balancer_subgrp(X, trt, Z, weightfunc, weightptr,
                                proxfunc, hyperparam, ridge, Q, NULL, NULL, opts)
    } else if(type == "subgrp_multi") {
        out <- balancer_multi(X, V, trt, Z, weightfunc, weightptr,
                              proxfunc, hyperparam, ridge, interact, opts)
    } else if(type == "multimatch") {
        out <- balancer_multimatch(X, V, trt, Z, weightfunc, weightptr,
                              proxfunc, hyperparam, ridge, opts)
    } else if(type == "missing") {
        out <- balancer_missing(X, trt, Z, weightfunc, weightptr,
                                proxfunc, hyperparam, ridge, Q, opts)
    } else if(type == "hte") {
        opts$linear <- linear
        out <- balancer_hte(X, trt,weightfunc, weightptr,
                            proxfunc, hyperparam, ridge, Q, opts)
    } else {
        stop("type must be one of ('att', 'subgrp', 'missing', 'hte')")
    }

    return(out)
}



balancer_subgrp <- function(X, trt, Z=NULL, weightfunc, weightfunc_ptr,
                             proxfunc, hyperparam, ridge, Q=NULL, kernel, kern_param=1, opts=list()) {
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
                     ridge=ridge
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
    for(i in 1:m) {
        k = grps[i]
        weights[(trt == 0) & (Z == k)] <-
            weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE], theta[,i,drop=FALSE])
    }

    out$weights <- weights

    ## The final imbalance
    loss_opts$ridge <- F
    out$imbalance <- balancing_grad(theta, loss_opts)
    return(out)
}



balancer_missing <- function(X, trt, R, weightfunc, weightfunc_ptr,
                             proxfunc, hyperparam, ridge, Q, opts=list()) {
    #' Balancing weights for missing outcomes
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param R vector of missingness indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param hyperparam Regularization hyper parameter
    #' @param ridge Whether to use L2 penalty in dual
    #' @param Q m x m matrix to tie together ridge penalty, default: NULL,
    #'          if TRUE, use covariance of treated groups
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
    x_t <- cbind(colMeans(X[trt==1,, drop=FALSE]),
                 colMeans(X[trt==1,, drop=FALSE]))

    obs_trt <- trt[R==1]
    loss_opts = list(Xc=X[R==1,,drop=FALSE],
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     weight_type="missing",
                     trt=obs_trt,
                     ridge=ridge
                     )

    loss_opts$idx_ctrl = which(obs_trt == 0) - 1
    loss_opts$idx_trt = which(obs_trt == 1) - 1

    ## if ridge, make the identity matrix
    if(ridge) {
        loss_opts$hyper <- hyperparam
        if(!is.null(Q)) { 
            loss_opts$hasQ <- TRUE
            if(Q) {
                loss_opts$Q <- t(x_t)
            } else {
                loass_opts$Q <- Q
            }
        } else {
            loss_opts$hasQ <- FALSE
        }
        
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
    init = matrix(0, nrow=dim(X)[2], ncol=m)
    
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))
    

    apgout <- apg(make_balancing_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel, opts$verbose)
                  

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




balancer_hte <- function(X, trt, weightfunc, weightfunc_ptr,
                          proxfunc, hyperparam, ridge, Q, opts=list()) {
    #' Balancing weights for heterogeneous treatment effects
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param hyperparam Regularization hyper parameter
    #' @param ridge Whether to use L2 penalty in dual
    #' @param Q m x m matrix to tie together ridge penalty, default: NULL,
    #'          if TRUE, use covariance of treated groups
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
    x_t <- t(X[trt==1,, drop=FALSE])

    ## initialize at 0
    init = matrix(0, nrow=dim(X)[2], ncol=m)

    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))
    
    loss_opts = list(Xc=X[trt==0,,drop=FALSE],
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     weight_type="hte",
                     linear=opts$linear,
                     ridge=ridge
                     )

    ## if ridge, make the identity matrix
    if(ridge) {
        loss_opts$hyper <- hyperparam
        if(!is.null(Q)) { 
            loss_opts$hasQ <- TRUE
            if(Q) {
                loss_opts$Q <- X[trt==1,,drop=FALSE]
            } else {
                loass_opts$Q <- Q
            }
        } else {
            loss_opts$hasQ <- FALSE
        }
    }

    
    if(length(hyperparam) == 1) {
        prox_opts = list(lam=hyperparam)
    } else if(length(hyperparam) == 2) {
        prox_opts = list(lam1=list(lam=hyperparam[1]),
                         lam2=list(lam=hyperparam[2]))
    } else {
        stop("hyperparam must be length 1 or 2")
    }

    apgout <- apg(make_balancing_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel, opts$verbose)
                  

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
