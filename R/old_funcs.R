################################################################################
## Dumping ground for old and unused stuff, probably get rid of this...
################################################################################



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

