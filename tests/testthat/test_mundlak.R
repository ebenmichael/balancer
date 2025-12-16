context("Testing that mundlak weights run")

set.seed(1011)
sample_simplex <- function(n) {

    u <- c(0, sort(runif(n - 1)), 1)
    return(diff(u))

}



## make fake data
n <- 1000
d <- 2
k <- 10
X <- matrix(rnorm(n * d), nrow = n)
z_probs <- sample_simplex(k)
Z <- sample(1:k, n, replace = T, prob = z_probs)
k <- length(unique(Z))
pscore <- 1 / (1 + exp(-rowSums(X)))
trt <- sapply(1:n, 
  function(i) {
    sample(c(0, 1), 1, prob = c(1 - pscore[i], pscore[i]))
  })

ind_covs <- X
clus_covs <- model.matrix(~ as.factor(Z) - 1) %*% t(sapply(split.data.frame(cbind(X, trt, X * trt), as.factor(Z)), colMeans))
interact_covs <- model.matrix(~ ind_covs : clus_covs - 1)

wts <- mundlak_weights(ind_covs, clus_covs, interact_covs, trt, Z, lambda = 1e10, 0, 1, T, T, T, F)
