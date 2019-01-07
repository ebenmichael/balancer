################################################################################
## Collection of balance functions corresponding to Regularization
################################################################################

#' Infinity norm
linf <- function(x) norm(matrix(x), "I") 


#' Two norm
l2 <- function(x) norm(matrix(x), "2")

#' Two norm squared
l2sq <- function(x) norm(matrix(x), "2")^2/2

#' Q norm
l2Q <- function(x, Q) {as.numeric(sqrt( t(x) %*% Q %*% x))}


#' One norm
l1 <- function(x) norm(matrix(x), "O") 

#' Group infinity norm
linf_grp <- function(x) sum(apply(x, 2, function(y) max(abs(y))))

#' Group inf + inf
linf_grp_linf <- function(x) {
    d <- ncol(x)
    linf(x[,1:d]) + linf_grp(x[,(d+1):ncol(x)])
}

#' Operator norm
op_norm <- function(x) norm(x, "2")

#' Operator + inf
linf_op <- function(x) {
    d <- ncol(x)
    linf(x[,1:d]) + op_norm(x[,(d+1):ncol(x)])
}
