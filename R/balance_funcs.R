################################################################################
## Collection of balance functions corresponding to Regularization
################################################################################

#' Infinity norm
linf <- function(x) norm(matrix(x), "I") 


#' Two norm
l2 <- function(x) norm(matrix(x), "2")

#' One norm
l1 <- function(x) norm(matrix(x), "O") 
