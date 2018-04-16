################################################################################
## Functions to solve the dual problem
################################################################################

#### Helper weights functions
softmax <- function(eta) {
    #' Compute numerically stable softmax with natural param eta
    #' Logit link
    m <- max(eta)
    num <- as.numeric(exp(eta - m))
    denom <- sum(exp(eta - m))
    return(num / denom)

}


normlin <- function(eta) {
    #' Linear odds, normalized weights
    return(eta / (sum(eta) + 0.00001))
    
}


poslin <- function(eta) {
    #' Linear odds, constrained to be positive
    return(eta * (eta >= 0))
}


##### Helper prox functions
no_prox <- function(x, lam) {
    #' Prox of 0 is the identity
    return(x)
}


l1_prox <- function(x, lam) {
    #' L1 norm prox operator, soft thresholding
    #' @param x input
    #' @param lam scaling function
    #'
    #' @return result of prox
    out <- (x - lam) * (x > lam) + (x + lam) * (x < -lam)
    return(out)
    
}


linf_prox <- function(x, lam) {
    #' L infinity norm prox operator
    #' Uses algorithm for projection onto the L1 ball from Duchi (2008)
    #' @param x input
    #' @param lam scaling function
    #'
    #' @return result of prox

    ## compute projection onto L1 ball
    absx <- abs(x)
    if(sum(absx) <= lam) {
        w <- x
    } else {
        absx <- sort(absx, decreasing = TRUE)
        
        ## get cumulative mean
        csx <- cumsum(absx)
        
        rho <- max(which(absx - (csx - lam) /  1:length(x) > 0))
        theta <- (csx[rho]- lam) / rho

        w <- ((x - theta) > 0)  * (x - theta)
        w <- sign(x) * w
    }
    return(x-w)
}


l1_grp_prox <- function(x, lam) {
    #' Group L1 norm prox operator, soft thresholding
    #' @param x input
    #' @param lam scaling function
    #'
    #' @return result of prox

    ## helper shrink function
    shrinkfunc <- function(xrow, lam) {
        max(0, 1 - lam / norm(xrow, "2"))
    }

    ## soft threshold along groups
    t(apply(x, 1, function(xrow) shrinkfunc(xrow, lam) * xrow))
}


nuc_prox <- function(x, lam) {
    #' Nuclear norm prox operator, soft thresholding
    #' @param x input
    #' @param lam scaling function
    #'
    #' @return result of prox

    ## compute svd
    ##out <- RSpectra::svds(x, 10)#min(dim(x))/2)
    out <- svd(x)
    ## soft threshold singular values
    s <- l1_prox(out$d, lam)

    ## combine back into matrix
    return(out$u %*% diag(s) %*% t(out$v))
}


balancer_subgrp <- function(X, trt, Z=NULL, weightfunc, proxfunc, hyperparam,
                         normalized=TRUE, opts=list()) {
    #' Helper function to fit the dual for general odds function and prox
    #' estimates heterogeneous treatment effects
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators
    #' @param weightfunc Derivative of convex conjugate of dispersion function (possibly normalized)
    #' @param proxfunc Prox operator of regularization function
    #' @param hyperparam Regularization hyper parameter
    #' @param normalized Whether to fit normalized weights, default: True
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
    grad <- function(theta, ...) {
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

    prox <- function(theta, step, ...) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)
        if(normalized) {
            ## apply prox operator for covariate parameters
            proxtheta <- proxfunc(theta[-1,,drop=FALSE], step * hyperparam)
            proxtheta <- matrix(proxtheta, ncol=m)
            ## prox is identity for bias
            proxtheta <- rbind(theta[1,,drop=FALSE], proxtheta)
        } else {
               ## apply prox operator for covariate parameters
            proxtheta <- proxfunc(theta, step * hyperparam)
        }
        
        return(as.numeric(proxtheta))
    }

    ## combine opts with defauls
    opts <- c(opts, MAX_ITERS=5000, EPS=1e-8)
    apgout <- apg::apg(grad, prox, d * m, opts=opts)

    ## collect results
    out <- list()

    ## theta
    theta <- matrix(apgout$x, ncol=m)
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


balancer_missing <- function(X, trt, R, weightfunc, proxfunc,
                             hyperparam, opts=list()) {
    #' Helper function to fit the dual for general odds function and prox
    #' Estimates ATT with missing outcomes
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param R vector of missingness indicators
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

    ## one set of weights for R=1,T=1 -> T=1
    ## one for R=1, T=0 -> T=1
    m <- 2

    n <- dim(X)[1]

    ## if(normalized) {
    ##     ## add a bias term
    ##     X <- cbind(rep(1, n), X)
    ## }

    d <- dim(X)[2]

    ## get the moments for treated units twice
    x_t <- cbind(colMeans(X[trt==1,]),
                 colMeans(X[trt==1,]))


    ## helper function to get weights for a given theta
    weights_func <- function(theta) {
        
        ## weights for R=1 T=0 to T=1
        weights1 <- numeric(n)
        weights1[trt == 0 & R == 1] <-
            weightfunc(X[trt == 0 & R == 1,,drop=FALSE] %*% theta[,1])
        
        ## weights1 <- sapply(1:n, function(i) {
        ##     if(trt[i] == 0 & R[i] == 1) {
        ##         conjprime(t(X[i,]) %*% theta[,1])
        ##     } else {
        ##         0
        ##     }})
        
        ## weights for R=1 T=1 to T=1
        weights2 <- numeric(n)
        weights2[trt == 1 & R == 1] <-
            weightfunc(X[trt == 1 & R == 1,,drop=FALSE] %*% theta[,1])

        ## weights2 <- sapply(1:n, function(i) {
        ##     if(trt[i] == 1 & R[i] == 1) {
        ##         conjprime(t(X[i,]) %*% theta[,2])
        ##     } else {
        ##         0
        ##     }})
       
        cbind(weights1, weights2)
    }

    
    ## objective gradient
    ## f^*'(x * theta) - x
    grad <- function(theta, ...) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)

        ## first part of gradient comes from units with observed outcomes
        ## for trt and ctrl
        weights <- weights_func(theta)

        grad1 <- t(X) %*% weights
        
        ## second part of gradient comes from all ctrl and treated units
        grad <- grad1 - x_t

        return(as.numeric(grad))
    }

    ## prox operator of regularizer

    prox <- function(theta, step, ...) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)

        ## if(normalized) {
        ##     ## apply prox operator for covariate parameters
        ##     proxtheta <- proxfunc(theta[-1,], step * hyperparam)
            
        ##     ## prox is identity for bias
        ##     proxtheta <- rbind(theta[1,], proxtheta)
        ## } else {
               ## apply prox operator for covariate parameters
            proxtheta <- proxfunc(theta, step * hyperparam)
        ## }
        
        return(as.numeric(proxtheta))
    }

    ## combine opts with defauls
    opts <- c(opts, MAX_ITERS=5000, EPS=1e-8)
    apgout <- apg::apg(grad, prox, d * m, opts=opts)

    ## collect results
    out <- list()

    ## theta
    theta <- matrix(apgout$x, ncol=m)

    out$theta <- theta

    out$weights <- weights_func(theta)

    ## The final imbalance
    out$imbalance <- t(X) %*% out$weights - x_t

    return(out)
}



balancer_individ <- function(X, trt, weightfunc, proxfunc,
                             hyperparam, opts=list()) {
    #' Helper function to fit the dual for general odds function and prox
    #' Estimates a CATE for each treated individual
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
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


    ## one set of weights for each treated unit
    m <- sum(trt)

    n <- dim(X)[1]

    ## if(normalized) {
    ##     ## add a bias term
    ##     X <- cbind(rep(1, n), X)
    ## }

    d <- dim(X)[2]
    
    ## keep the covariates for the treated units
    x_t <- t(X[trt==1,])

    ## covariates for control units
    x_c <- X[trt==0,]

    ## helper function to get weights for a given theta
    weights_func <- function(theta) {
        eta <- X %*% theta
        eta[trt==1,] <- 0
        weights <- apply(eta, 2, weightfunc)
        ## zero out weights on treated units
        ## weights[trt==1,] <- 0
        weights
    }

    
    ## objective gradient
    ## f^*'(x * theta) - x
    grad <- function(theta, ...) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)

        ## first part of gradient comes from units with observed outcomes
        ## for trt and ctrl
        weights <- weights_func(theta)

        grad1 <- t(X) %*% weights

        ## second part of gradient comes from all ctrl and treated units
        grad <- grad1 - x_t
        return(as.numeric(grad))
    }

    ## prox operator of regularizer

    prox <- function(theta, step, ...) {
        ## reshape paramters into matrix
        theta <- matrix(theta, ncol=m)

        ## if(normalized) {
        ##     ## apply prox operator for covariate parameters
        ##     proxtheta <- proxfunc(theta[-1,], step * hyperparam)
            
        ##     ## prox is identity for bias
        ##     proxtheta <- rbind(theta[1,], proxtheta)
        ## } else {
               ## apply prox operator for covariate parameters
        proxtheta <- proxfunc(theta, step * hyperparam)
        ## }
        
        return(as.numeric(proxtheta))
    }

    ## combine opts with defauls
    opts <- c(opts, MAX_ITERS=5000, EPS=1e-8)
    apgout <- apg::apg(grad, prox, d * m, opts=opts)

    ## collect results
    out <- list()

    ## theta
    theta <- matrix(apgout$x, ncol=m)

    out$theta <- theta

    out$weights <- weights_func(theta)

    ## The final imbalance
    out$imbalance <- t(X) %*% out$weights - x_t

    return(out)
}
