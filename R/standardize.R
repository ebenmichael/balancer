################################################################################
## Wrapper to standardize to target means
################################################################################




standardize <- function(X, target,
                        Z=NULL, 
                        link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                        regularizer=c("ridge", "l1", "grpl1", "l2", "linf", "nuc",
                                      "l1_all", "l1_nuc"),
                        lambda=NULL, nlambda=20, lambda.min.ratio=1e-3,
                        normalized=TRUE, alpha=1,
                        opts=list()) {
    #' Re-weight groups to target population means
    #' @param X n x d matrix of covariates
    #' @param target Vector of population moments to weight towards
    #' @param Z Vector of site indicators
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param lambda Regularization hyperparameter
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param interact Whether to interact group and individual level covariates
    #' @param normalized Whether to normalize the weights
    #' @param alpha Elastic net parameter \eqn{\frac{1-\alpha}{2}\|\beta\|_2^2 + \alpha\|\beta\|_1}, defaults to 1
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
    params <- map_to_param(X, link, regularizer, NULL, normalized, Q, alpha)
    weightfunc <- params[[1]]
    weightptr <- params[[2]]
    proxfunc <- params[[3]]
    balancefunc <- params[[4]]
    prox_opts <- params[[5]]


    prep <- preprocess(X, trt, NULL, "att", link, normalized)
    X <- prep$X
    ipw_weights <- prep$ipw_weights

    if(normalized) target <- c(1, target)

    out <- standardize_(X, target, Z, weightfunc, weightptr, proxfunc,
                        balancefunc, lambda, nlambda, lambda.min.ratio,
                        ipw_weights, opts, prox_opts)
    return(out)
}


standardize_  <- function(X, target, Z, weightfunc, weightfunc_ptr,
                            proxfunc, balancefunc, lambda, nlambda,
                            lambda.min.ratio, ipw_weights, 
                          opts=list(), prox_opts) {
    #' Internal function to reweight to target population means
    #' @param X n x d matrix of covariates
    #' @param target Vector of population means to re-weight to
    #' @param Z Vector of hierarchical factor indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param balancefunc Balance criterion measure
    #' @param lambda Regularization hyper parameter
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}}
    #' @param prox_opts List of additional arguments for prox    
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}
    
    ## if no subgroups, put everything into one group
    if(is.null(Z)) {
        Z <- rep(0, length(trt))
    }


    
    ## get the distinct group labels
    grps <- sort(unique(Z))
    m <- length(grps)

    n <- dim(X)[1]
    d <- dim(X)[2]
    
    target <- as.matrix(sapply(grps,
                  function(k) target))


    loss_opts = list(X=X,
                     target=target,
                     z=Z,
                     weight_func=weightfunc_ptr,
                     ipw_weights=ipw_weights,
                     n_groups=m,
                     dim=d
                     )


    ## initialize at 0
    init = matrix(0, nrow=d, ncol=(m))
    
    ## combine opts with defauls
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))

    ## if hyperparam is NULL, start from reference weights and decrease
    if(is.null(lambda)) {
        lam0 <- balancefunc(balancing_grad_multilevel(init, loss_opts))
        lam1 <- lam0 * lambda.min.ratio
        ## decrease on log scale
        lambda <- exp(seq(log(lam0), log(lam1), length.out=nlambda))
    }


    ## collect results
    out <- list()
    out$theta <- matrix(,nrow=d, ncol=length(lambda))
    out$imbalance <- matrix(,nrow=d, ncol=length(lambda))    
    out$weights <- matrix(0, nrow=n, ncol=length(lambda))
    out$weightfunc <- weightfunc


    ## with multiple hyperparameters do warm starts        
    prox_opts <- c(prox_opts,
                  list(lam=1))


    apgout <- apg_warmstart(make_balancing_grad_standardize(),
                            proxfunc, loss_opts, prox_opts,
                            lambda,
                            opts$x, opts$max_it, opts$eps,
                            opts$alpha, opts$beta, opts$accel, opts$verbose)

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    weights_from_theta <- function(th) {
        weights <- matrix(0, nrow=n, ncol=m)
        for(i in 1:m) {
            k = grps[i]
            weights[(Z == k), i] <-
                weightfunc(X[Z == k, , drop=FALSE],
                           th[,i,drop=FALSE],
                           ipw_weights[Z == k,,drop=FALSE])
        }
        weights
    }
    
    weights <- lapply(out$theta, weights_from_theta)

    out$weights <- if(length(weights)==1) weights[[1]] else weights

    
    ## The final imbalance    
    imbalance <- lapply(theta, function(th) balancing_grad_standardize(th, loss_opts))
    out$imbalance <- if(length(imbalance) == 1) imbalance[[1]] else imbalance

    out$theta <- if(length(out$theta) == 1) out$theta[[1]] else out$theta

    out$lambda <- lambda
    return(out)
}
