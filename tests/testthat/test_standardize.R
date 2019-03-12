context("Testing that standardization runs")

## make fake data
n <- 1000
d <- 100
k <- 20
X <- matrix(rnorm(n * d), nrow=n)
Z <- sample(1:k, n, replace=T)
target <- colMeans(X)
                            
test_that("Standardization runs without hiccups", {
    
    out <- standardize(X, target, Z, lambda=1e1)

    ## weights are right shape
    expect_equal(dim(out$weights), c(n, k))

    ## dual parameters are right shape
    expect_equal(dim(out$theta), c(d+1, k))
    
    
}
)
