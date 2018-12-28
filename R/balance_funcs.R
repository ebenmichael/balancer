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
l2Q <- function(x, Q) {sqrt( t(x) %*% Q %*% x)}


#' One norm
l1 <- function(x) norm(matrix(x), "O") 
