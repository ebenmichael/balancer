context("Testing that standardization runs")

## make fake data
n <- 1000
d <- 10
k <- 20
X <- matrix(rnorm(n * d), nrow=n)
Z <- sample(1:k, n, replace=T)
target <- colMeans(X)
                            
test_that("Standardization runs without hiccups", {
    
    out <- standardize(X, target, Z, lambda = 0, verbose = F)

    ## weights are right shape
    expect_equal(dim(out$weights), c(n, k))

    ## imbalance is the right shape
    expect_equal(dim(out$imbalance), c(d, k))

    ## P matrix is the right shape
    expect_equal(dim(out$data_out$P), c(n, n))
}
)
