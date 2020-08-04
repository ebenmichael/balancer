library(dplyr)
source("../R/standardize_treatment.R")
source("../R/standardize_treatment_hybrid.R")
source("../R/standardize_treatment_kernel.R")

################################################################################
### generate data and optimization inputs
################################################################################

# sample size in each site
n1 <- 10
n2 <- 16
n3 <- 20
n <- n1 + n2 + n3

# site indicator
S <- c(rep(1, n1), rep(2, n2), rep(3, n3))

# covariates
X <- data.frame(V1 = rnorm(n, mean = S, sd = S),
                V2 = rbinom(n, 1, (4 - S) / 4))
X$V1[S == 1] <- c(1.2, 1.2, 1.2, 1.2, -1, -1, 1.2, 1.2, 1.2, 1.2)
X$V2[S == 1] <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0)

# propensity scores
pscores <- c(0.5, 0.5, 0.75)

# treatment indicators
Z <- c(sample(c(rep(1, pscores[1] * n1), rep(0, (1 - pscores[1]) * n1))),
       sample(c(rep(1, pscores[2] * n2), rep(0, (1 - pscores[2]) * n2))),
       sample(c(rep(1, pscores[3] * n3), rep(0, (1 - pscores[3]) * n3))))

Z[1:n1] <- rep(c(1, 0), n1 / 2)

pscores <- rep(pscores, c(n1, n2, n3))

df <- data.frame(S, X, pscores, Z)
df2 <- df #sample_frac(df)
S <- df2$S
X <- df2[, 2:3]
pscores <- df2$pscores
Z <- df2$Z

# target
target <- c(0.1)

# arguments
lambda <- 0
lowlim <- 0
uplim <- max(c(n1, n2, n3))
scale_sample_size <- T
data_in <- NULL
verbose <- TRUE
return_program <- TRUE
exact_global <- F
init_uniform <- F
eps_abs <- 1e-5
eps_rel <- 1e-5
kernel0 <- kernlab::vanilladot()
kerneltau <- kernlab::vanilladot()
gc <- TRUE

X0 <- X
Xtau <- X[, 1]

# remove excess variables
rm(n)
rm(X)

################################################################################
### compare outputs
################################################################################

a <- standardize_treatment_hybrid(X0 = X0,
                                  Xtau = Xtau,
                                  target = target,
                                  S = S,
                                  Z = Z,
                                  pscores = pscores,
                                  kernel0 = kernlab::vanilladot(),
                                  lambda = lambda,
                                  lowlim = lowlim,
                                  uplim = uplim,
                                  scale_sample_size = scale_sample_size,
                                  data_in = NULL,
                                  verbose = TRUE,
                                  return_program = FALSE,
                                  gc = TRUE)

b <- standardize_treatment(X0 = X0,
                           Xtau = Xtau,
                           target = target,
                           S = S,
                           Z = Z,
                           pscores = pscores,
                           lambda = lambda,
                           lowlim = lowlim,
                           uplim = uplim,
                           scale_sample_size = scale_sample_size,
                           data_in = NULL,
                           verbose = TRUE,
                           return_program = FALSE)

# they're not the same
a$weights
b$weights

################################################################################
### preprocessing in both hybrid and linear versions
################################################################################

# ensure that covariate matrices are matrices and get total number of units
X0 <- as.matrix(X0)
Xtau <- as.matrix(Xtau)
n <- nrow(X0)

# convert site indicators to factor, get index of sites, and compute # sites
S_factor <- as.factor(S)
unique_S <- levels(S_factor)
J <- length(unique_S)

# split covariate matrix by site and compute number of units in each site
X0s <- split.data.frame(X0, S_factor)
Xtaus <- split.data.frame(Xtau, S_factor)
nj <- sapply(X0s, nrow)

# get row IDs of units by site
idxs <- split(1:nrow(X0), S_factor)

# ensure that target is a vector
target <- c(target)

# create propensity score multipliers and split by site
pro_trt <- Z / (pscores * nj[S_factor])
pro_trt_split <- split(pro_trt, S_factor)
pro_ctr <- (1 - Z) / ((1 - pscores) * nj[S_factor])
pro_ctr_split <- split(pro_ctr, S_factor)

# define number of auxiliary weights (only the linear version uses this)
aux_dim <- (J * ncol(X0)) + (J * ncol(Xtau))

################################################################################
### create q vector
################################################################################

# linear
q_linear <- create_q_vector_treatment(Xtaus, pro_trt_split, target, aux_dim)
# hybrid
q_hybrid <- create_q_vector_treatment(Xtaus, pro_trt_split, target, 0)

################################################################################
### create p matrix
################################################################################

# linear
P_linear <- create_P_matrix_treatment(n, aux_dim)
# hybrid
P_hybrid <- create_P_matrix_treatment_kernel(n, X0s, Xtaus, kernel0, kernlab::vanilladot(),
                                             pro_trt, pro_ctr, S_factor, gc)

################################################################################
### create arbitrary weights (to check that objective function is the same)
################################################################################

gamma1 <- runif(n1, 0, 1)
gamma1 <- n1 * gamma1 / sum(gamma1)
gamma2 <- runif(n2, 0, 1)
gamma2 <- n2 * gamma2 / sum(gamma2)
gamma3 <- runif(n3, 0, 1)
gamma3 <- n3 * gamma3 / sum(gamma3)
gamma <- c(gamma1, gamma2, gamma3)

J <- length(X0s)
nj <- as.numeric(sapply(X0s, nrow))
d0 <- ncol(X0s[[1]])
dtau <- ncol(Xtaus[[1]])
n <- sum(nj)
X0st <- lapply(X0s, t)
Xtaust <- lapply(Xtaus, t)

# compute number of treated and control units within each site
n1j <- sapply(unique(S_factor), FUN = function(j) sum(Z[S_factor == j]))
n0j <- sapply(unique(S_factor), FUN = function(j) sum(1 - (Z[S_factor == j])))

sqrtP1 <- Matrix::t(Matrix::t(Matrix::bdiag(X0st)) * (unlist(pro_trt_split) - unlist(pro_ctr_split)))
sqrtP2 <- Matrix::t(Matrix::t(Matrix::bdiag(Xtaust)) * unlist(pro_trt_split))

gamma_aux <- c(as.vector(sqrtP1 %*% gamma), as.vector(sqrtP2 %*% gamma))

################################################################################
### check that objective function is the same
################################################################################

# linear
0.5 * (c(gamma, gamma_aux) %*% P_linear %*% c(gamma, gamma_aux)) + q_linear %*% c(gamma, gamma_aux)
# hybrid
0.5 * (gamma %*% P_hybrid %*% gamma) + q_hybrid %*% gamma
