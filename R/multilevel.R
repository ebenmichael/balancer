################################################################################
## Balancer for multilevel modelling
################################################################################



balancer_multi  <- function(X, V=NULL, trt, Z, weightfunc, weightfunc_ptr, proxfunc,
                            hyperparam, ridge=F, opts=list()) {
    #' Balancing weights for ATT (with hierarchical structure)
    #' @param X n x d matrix of covariates
    #' @param V n x k matrix of group level covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of hierarchical factor indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
    #' @param hyperparam Regularization hyper parameter
    #' @param ridge Whether to use ridge penalty in dual
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

    if(ridge) {
        global_param = hyperparam[3]
        group_param = hyperparam[4]
    } else {
        global_param = 0
        group_param = 0
    }
    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     n_groups=m,
                     dim=d,
                     global_int=hyperparam[1],
                     group_int=hyperparam[2],
                     global_param=global_param,
                     group_param=group_param,
                     weight_type="subgroup",
                     z=Z[trt==0],
                     z_ind = lapply(grps, function(k) which(ctrl_z==k)-1),
                     grp_cov_dim=grp_cov_dim,
                     ridge=ridge
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

    if(length(hyperparam) == 1) {
        prox_opts = list(lam=hyperparam)
    } else if(length(hyperparam) == 4) {
        prox_opts = list(lam1=list(lam=hyperparam[1]),
                         lam2=list(lam=hyperparam[2]),
                         lam3=list(lam=hyperparam[3]),
                         lam4=list(lam=hyperparam[4]),
                         n_groups=m)
    } else {
        stop("hyperparam must be length 1 or 4")
    }
    
    apgout <- apg(make_multilevel_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel, opts$verbose)

    ## collect results
    out <- list()

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    weights <- matrix(0, nrow=n, ncol=m)
    for(i in 1:m) {
        k = grps[i]
        weights[(trt == 0) & (Z == k), i] <-
            weightfunc(X[(trt == 0) & (Z == k), , drop=FALSE], theta[,1,drop=FALSE] + theta[,i+1,drop=FALSE])
    }
    

    out$weights <- weights

    ## The final imbalance

    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     n_groups=m,
                     dim=d,
                     global_int=0,
                     group_int=0,
                     global_param=0,
                     group_param=0,
                     weight_type="subgroup",
                     z=Z[trt==0],
                     z_ind = lapply(grps, function(k) which(ctrl_z==k)-1),
                     grp_cov_dim=grp_cov_dim,
                     ridge=ridge
                     )
    
    out$imbalance <- multilevel_grad(theta, loss_opts)
    return(out)
}


balancer_multimatch  <- function(X, V=NULL, trt, Z, weightfunc, weightfunc_ptr, proxfunc,
                            hyperparam, opts=list()) {
    #' Balancing weights for ATT with multilevel treatment assignment
    #' @param X n x d matrix of covariates
    #' @param V n x k matrix of group level covariates    
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of hierarchical group indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param weightfunc_ptr Pointer to weightfunc
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


    ## get treated group labels
    grps <- sort(unique(Z[trt==1]))
    m <- length(grps)
    
    n <- dim(X)[1]
    d <- dim(X)[2]

    X <- cbind(1, X)


    ## append the group level covariates
    if(!is.null(V)) {
        X <- cbind(X, V)
    }
    
    ## get treated means in treated groups
    x_t <- sapply(grps,
                  function(k) colSums(X[(trt ==1) & (Z==k), , drop=FALSE]))
    x_t <- as.matrix(x_t)
    x_t <- cbind(colSums(X[trt==1,]), x_t) / sum(trt==1)
        
    Xc <- X[trt==0,,drop=FALSE]


    grp_cov_dim <- if(is.null(V)) 0 else dim(V)[2]

    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     n_groups=m,
                     dim=d,
                     global_int=hyperparam[1],
                     group_int=hyperparam[2],
                     global_param=hyperparam[3],
                     group_param=hyperparam[4],
                     weight_type="multimatch",
                     grp_cov_dim=grp_cov_dim
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


    if(length(hyperparam) == 1) {
        prox_opts = list(lam=hyperparam)
    } else if(length(hyperparam) == 4) {
        prox_opts = list(lam1=list(lam=hyperparam[1]),
                         lam2=list(lam=hyperparam[2]),
                         lam3=list(lam=hyperparam[3]),
                         lam4=list(lam=hyperparam[4]),
                         n_groups=m)
    } else {
        stop("hyperparam must be length 1 or 4")
    }
    
    apgout <- apg(make_multilevel_grad(), proxfunc, loss_opts, prox_opts,
                  opts$x, opts$max_it, opts$eps, opts$alpha, opts$beta, opts$accel, opts$verbose)
    ## collect results
    out <- list()

    ## theta
    theta <- apgout
    out$theta <- theta
    ## weights
    weights <- matrix(0, nrow=n, ncol=m)
    for(i in 1:m) {
        k = grps[i]
        weights[(trt == 0), i] <-
            weightfunc(X[(trt == 0), , drop=FALSE], theta[,1,drop=FALSE] + theta[,i+1,drop=FALSE])
    }
    

    out$weights <- weights

    ## The final imbalance

    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     weight_func=weightfunc_ptr,
                     n_groups=m,
                     dim=d,
                     global_int=0,
                     group_int=0,
                     global_param=0,
                     group_param=0,
                     weight_type="multimatch",
                     grp_cov_dim=grp_cov_dim
                     )
    
    out$imbalance <- multilevel_grad(theta, loss_opts)
    return(out)
}
