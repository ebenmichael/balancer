context("Testing that multilevel_qp runs")

set.seed(1011)
sample_simplex <- function(n) {

    u <- c(0, sort(runif(n - 1)), 1)
    return(diff(u))

}



## make fake data
n <- 2500
d <- 1
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
X_fixed_eff <- model.matrix(~ as.factor(Z) + X - 1)


test_that("Two different ways of ignoring local balance are equivalent", {
  out1 <- multilevel_qp(X_fixed_eff, trt, rep(1, length(trt)),
                         lambda = 1e10, verbose = F)
  out2 <- multilevel_qp(X, trt, Z,
                         lambda = 1e8, verbose = T, scale_sample_size = F)
  expect_equal(out1$weights, out2$weights, tol=1e-2)
})
