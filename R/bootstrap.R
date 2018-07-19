################################################################################
## Bootstrap confidence intervals
################################################################################


balancer_boot <- function(nboot, X, trt, Z=NULL, type=c("att", "subgrp", "missing", "hte"),
                          link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                          regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc"),
                          hyperparam, Q=NULL, kernel=NULL, kern_param=1, normalized=TRUE,
                          sample_trt=FALSE, opts=list()) {
    #' Get point estimate and bootstraped standard errors
    #' @param nboot Number of bootstrap samples
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param Z Vector of subgroup indicators or observed indicators
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             ATT with missing outcomes, and heterogeneouts effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param hyperparam Regularization hyperparameter
    #' @param Q m x m matrix to tie together ridge penalty, default: NULL,
    #'          if TRUE, use covariance of treated groups
    #' @param kernel What kernel to use, default NULL
    #' @param kern_param Hyperparameter for kernel
    #' @param normalized Whether to normalize the weights
    #' @param sample_trt Whether to resample treated units
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

    
}
