################################################################################
## Balancer for multilevel modelling
################################################################################



balancer_multi  <- function(X, V=NULL, trt, Z, weightfunc, weightfunc_ptr,
                            proxfunc, balancefunc, lambda, nlambda,
                            lambda.min.ratio, ipw_weights, 
                            interact=F, opts=list(), prox_opts) {
    #' Balancing weights for ATT (with hierarchical structure)
    #' @param X n x d matrix of covariates
    #' @param V n x k matrix of group level covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of hierarchical factor indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param proxfunc Prox operator of regularization function
    #' @param balancefunc Balance criterion measure
    #' @param lambda Regularization hyper parameter
    #' @param nlambda Number of hyperparameters to consider
    #' @param lambda.min.ratio Smallest value of hyperparam to consider, as proportion of smallest
    #'                         value that gives the reference weights
    #' @param ipw_weights Separately estimated IPW weights to measure dispersion against, default is NULL
    #' @param interact Whether to interact group and individual level characterisitics
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
    

    ## ## add indicators and interaction terms
    ## X_fac <- matrix(model.matrix(~as.factor(Z)-1, data.frame(X)),
    ##                 nrow=nrow(X))
    ## X_interact <- matrix(model.matrix(~.:as.factor(Z) -1, data.frame(X)),
    ##             nrow=nrow(X))    
    ## X <- cbind(1,X_fac, X, X_interact)

    n <- dim(X)[1]
    d <- dim(X)[2]
    
    X <- cbind(1, X)

    
    ## append the group level covariates
    if(!is.null(V)) {

        if(interact) {
            ## add in group and individual level interactions
            interacts <- matrix(model.matrix(~X:V), nrow=nrow(X))
            V <- cbind(V, interacts)
        }
        X <- cbind(X, V)
    }
    
    x_t <- sapply(grps,
                  function(k) colSums(X[(trt ==1) & (Z==k), , drop=FALSE]))
    x_t <- as.matrix(x_t)
    x_t <- cbind(colSums(X[trt==1,]), x_t)

    ## print(x_t)
    Xc <- X[trt==0,,drop=FALSE]


    grp_cov_dim <- if(is.null(V)) 0 else dim(V)[2]

    ## indices of subgroups
    ctrl_z = Z[trt==0]
    ipw_weights <- ipw_weights[trt==0,,drop=F]

    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     z_c=ctrl_z,
                     z_t=Z[trt==1],
                     weight_func=weightfunc_ptr,
                     ipw_weights=ipw_weights,
                     n_groups=m,
                     dim=d,
                     z=Z[trt==0],
                     z_ind = lapply(grps, function(k) which(ctrl_z==k)-1),
                     grp_cov_dim=grp_cov_dim,
                     weight_type="sugrp",
                     ridge=F
                     )


    ## initialize at 0
    init = matrix(0, nrow=(d+grp_cov_dim+1), ncol=(m+1))



    
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
    prox_opts = c(prox_opts,
                  list(lam=1))


    apgout <- apg_warmstart(make_balancing_grad_multilevel(),
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
            weights[(trt == 0) & (Z == k), i] <-
                weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE],
                           th[,1,drop=FALSE] + th[,i+1,drop=FALSE],
                           ipw_weights[Z[trt==0] == k,,drop=FALSE]) *
                sum(Z[trt==1] == k)
        }
        weights
    }
    
    out$weights <- lapply(out$theta, weights_from_theta)

    ## The final imbalance    
    out$imbalance <- lapply(theta, function(th) balancing_grad_multilevel(th, loss_opts))

    out$lambda <- lambda
    return(out)
}
