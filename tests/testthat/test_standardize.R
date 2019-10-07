context("Testing that standardization runs")

set.seed(1011)
sample_simplex <- function(n) {

    u <- c(0, sort(runif(n - 1)), 1)
    return(diff(u))

}



## make fake data
n <- 2500
d <- 10
k <- 100
X <- matrix(rnorm(n * d), nrow = n)
z_probs <- sample_simplex(k)
# labs <- c(1, 2, 4, 5, 39, 42, 43, 58, 66, 80, 85, 98, 103,
#           4750, 4810, 4820, 4890, 4910, 4920, 4960, 4980, 5000, 5010, 5020, 5040, 5100)
Z <- sample(1:k, n, replace = T, prob = z_probs)
target <- colMeans(X)


                            
test_that("Standardization throws errors for malformed data", {

    Znew <- Z
    Znew[1] <- NA
    expect_error(standardize(X, target, Znew, lambda = 0, verbose = F),
                 "Grouping vector Z contains NA values.")

    Xnew <- X
    Xnew[1,1] <- NA
    expect_error(standardize(Xnew, target, Z, lambda = 0, verbose = F),
                 "Covariate matrix X contains NA values.")

    targetnew <- target
    targetnew[1] <- NA
    expect_error(standardize(X, targetnew, Z, lambda = 0, verbose = F),
                 "Target vector contains NA values.")

    # targetnew <- target[-1]
    # expect_error(standardize(X, targetnew, Z, lambda = 0, verbose = F),
    #              paste0("Target dimension (", d - 1,
    #                     ") is not equal to data dimension (",
    #                      d, ")."))
}
)


test_that("Standardization runs without hiccups", {

    out <- standardize(X, target, Z, lambda = 0, verbose = F)

    ## weights are right shape
    expect_equal(dim(out$weights), c(n, k))

    ## imbalance is the right shape
    expect_equal(dim(out$imbalance), c(d, k))

    ## P matrix is the right shape
    expect_equal(dim(out$data_out$P), c(n + d * k, n + d * k))
}
)

test_that("Standardization runs with re-inputted data", {

    out1 <- standardize(X, target, Z, lambda = 0, verbose = F)

    out2 <- standardize(X, target, Z, lambda = 0, verbose = F, 
                        data_in = out1$data_out)

    ## weights are right shape
    expect_equal(out1$weights, out2$weights, tolerance = 1e-1)
}
)