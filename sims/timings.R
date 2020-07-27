# load required packages
library(balancer)
library(dplyr)
library(kernlab)
library(foreach)
library(doParallel)
library(doRNG)
library(tictoc)

################################################################################
### set parameters
################################################################################

# data-generating parameters (noise_sd will be refined)
noise_sd <- 1000

# standardization parameters
lambdas <- exp(seq(-8, 8, length.out = 6))
scale_samp <- FALSE
constrain_prop <- TRUE
exact_globe <- FALSE

################################################################################
### preprocessing
################################################################################

# load dataset
df <- read.csv("WW_datareqest_20170228.csv")

# get sample size, vector of unique offices, and vector of site-level covariate
frac <- 1
offices <- unique(df$office)
study <- ifelse(offices < 200, "A", ifelse(offices < 300, "B", "C"))

# remove redundant outcome column
df <- df[, -4]

# rename outcome column to denote that it is the potential outcome under control
names(df)[3] <- "y0"

# convert outcome to employment indicator
# df$y0 <- as.numeric(df$y0 > 0)

# define site-level covariate for each unit in the dataframe
df$studyA <- as.numeric(df$office < 200)
df$studyB <- as.numeric(df$office >= 200 & df$office < 300)
df$studyC <- as.numeric(df$office >= 300)

# define unit-level covariate columns
unit_covariate_cols <- 4:26

# define coefficient on site-level covariate to define heterogeneous treatment
# effects
beta <- c(2000, 6000, 0)

# define vector of coefficients on unit-level covariates to define heterogeneous
# treatment effects
theta <- c(
  # totearn pezero pelt2499 pet7499 pege7500 agele24 age2534 age3544 agege45 white
  0,        900,   0,       0,      0,       500,    0,      0,      0,      1002,
  # black hispanic natam asian oracenw childct chnum1 chnum2 chnum3 ychldlt6 hsged
  458,    820,    100,  1075, 1461,   0,      876,   1203,  1380,  1075,    0,
  # ahist12 applcnt
  0,        0)

# define covariate columns related to treatment
treatment_covariate_cols <- unit_covariate_cols[theta != 0]

# compute conditional average treatment effect for each unit
df$cate <- as.vector(data.matrix(df[, -(1:3)]) %*% matrix(c(theta, beta)))


# bootstrap the units within each site
bs_df <- df %>%
  group_by(office) %>%
  sample_frac(frac, replace = TRUE)

n <- nrow(bs_df)

# generate individual treatment effect for each unit
bs_df$ite <- bs_df$cate + rnorm(n, 0, noise_sd)

# define potential outcome under treatment
bs_df$y1 <- bs_df$y0 + bs_df$ite

# randomly assign units to treatment and control according to proportions from
# original data (recall that bs_df is already randomly ordered)
# bs_df$p <- unlist(sapply(unique(bs_df$office), FUN = function(x) sample(bs_df$p[bs_df$office == x])))
bs_df$p <- df$p

# set observed outcomes based on treatment assignment
bs_df$y <- ifelse(bs_df$p == 0, bs_df$y0, bs_df$y1)

# standardize unit-level covariates
bs_df[, unit_covariate_cols] <- scale(bs_df[, unit_covariate_cols])

# compute target covariate means (should all be zero)
target_x <- colMeans(bs_df[, unit_covariate_cols])

# compute propensity scores
pscores <- c(rep(mean(bs_df$p[bs_df$studyA == 1]), sum(bs_df$studyA)),
             rep(mean(bs_df$p[bs_df$studyB == 1]), sum(bs_df$studyB)),
             rep(mean(bs_df$p[bs_df$studyC == 1]), sum(bs_df$studyC)))

# set regularization parameter
lambda <- 0

############################################################################
### benchmark methods
############################################################################
# weight units in each site so that weighted means of unit-level covariates
# equal target
tic()
output <- standardize_rct(X = bs_df[, unit_covariate_cols],
                          trt = bs_df$p,
                          target = target_x,
                          Z = bs_df$office,
                          lambda = lambda,
                          scale_sample_size = scale_samp,
                          constrain_trt_prop = constrain_prop,
                          exact_global = exact_globe,
                          verbose = FALSE)
toc()

# weight units within each site so that each site is approximately balanced
# with the target site
tic()
trt_std <- standardize(X = bs_df[bs_df$p == 1, unit_covariate_cols],
                       target = target_x,
                       Z = bs_df$office[bs_df$p == 1],
                       lambda = lambda,
                       scale_sample_size = scale_samp,
                       exact_global = exact_globe,
                       verbose = FALSE)
ctr_std <- standardize(X = bs_df[bs_df$p == 0, unit_covariate_cols],
                       target = target_x,
                       Z = bs_df$office[bs_df$p == 0],
                       lambda = lambda,
                       scale_sample_size = scale_samp,
                       exact_global = exact_globe,
                       verbose = FALSE)
toc()

# try the non-kernel version of the treatment effect weighting estimator
tic()
blah <- standardize_treatment(X0 = bs_df[, unit_covariate_cols],
                              Xtau = bs_df[, treatment_covariate_cols],
                              target = target_x[(treatment_covariate_cols - 3)],
                              S = bs_df$office,
                              Z = bs_df$p,
                              pscores = pscores,
                              lambda = lambda,
                              lowlim = 0,
                              uplim = nrow(bs_df),
                              scale_sample_size = scale_samp,
                              data_in = NULL,
                              verbose = FALSE,
                              return_program = FALSE)
toc()

############################################################################
### kernel methods
############################################################################

# all sites at once
tic("total time")
hybrid_std <- standardize_treatment_hybrid(X0 = bs_df[, unit_covariate_cols],
                                           Xtau = bs_df[, treatment_covariate_cols],
                                           target = target_x[(treatment_covariate_cols - 3)],
                                           S = bs_df$office,
                                           Z = bs_df$p,
                                           pscores = pscores,
                                           kernel0 = kernlab::rbfdot(),
                                           lambda = lambda,
                                           lowlim = 0,
                                           uplim = max(table(bs_df$office)),
                                           scale_sample_size = scale_samp,
                                           data_in = NULL,
                                           verbose = FALSE,
                                           return_program = FALSE,
                                           gc = TRUE)
toc(log = TRUE)

# each site sequentially
tic("total time")
for(office in unique(bs_df$office)) {
  tic(paste0("office ", office))
  hybrid_std <- standardize_treatment_hybrid(X0 = bs_df[bs_df$office == office, unit_covariate_cols],
                                             Xtau = bs_df[bs_df$office == office, treatment_covariate_cols],
                                             target = target_x[(treatment_covariate_cols - 3)],
                                             S = rep(office, sum(bs_df$office == office)),
                                             Z = bs_df$p[bs_df$office == office],
                                             pscores = pscores[bs_df$office == office],
                                             kernel0 = kernlab::rbfdot(),
                                             lambda = lambda,
                                             lowlim = 0,
                                             uplim = sum(bs_df$office == office),
                                             scale_sample_size = scale_samp,
                                             data_in = NULL,
                                             verbose = FALSE,
                                             return_program = FALSE,
                                             gc = TRUE)
  toc(log = TRUE)
}
toc(log = TRUE)

# do one site ncores times in parallel
ncores <- detectCores()
tic("total time")
output <- mclapply(rep(120, ncores), FUN = function(office) {
  standardize_treatment_hybrid(X0 = bs_df[bs_df$office == office, unit_covariate_cols],
                               Xtau = bs_df[bs_df$office == office, treatment_covariate_cols],
                               target = target_x[(treatment_covariate_cols - 3)],
                               S = rep(office, sum(bs_df$office == office)),
                               Z = bs_df$p[bs_df$office == office],
                               pscores = pscores[bs_df$office == office],
                               kernel0 = kernlab::rbfdot(),
                               lambda = lambda,
                               lowlim = 0,
                               uplim = sum(bs_df$office == office),
                               scale_sample_size = scale_samp,
                               data_in = NULL,
                               verbose = FALSE,
                               return_program = FALSE,
                               gc = TRUE)},
  mc.cores = ncores # change this to see how the number of cores used affects the timings
)
toc(log = TRUE)
closeAllConnections()
rm(output)

# do one site ncores times sequentially
tic("total time")
for (office in rep(120, ncores)) {
  output <- standardize_treatment_hybrid(X0 = bs_df[bs_df$office == office, unit_covariate_cols],
                                       Xtau = bs_df[bs_df$office == office, treatment_covariate_cols],
                                       target = target_x[(treatment_covariate_cols - 3)],
                                       S = rep(office, sum(bs_df$office == office)),
                                       Z = bs_df$p[bs_df$office == office],
                                       pscores = pscores[bs_df$office == office],
                                       kernel0 = kernlab::rbfdot(),
                                       lambda = lambda,
                                       lowlim = 0,
                                       uplim = sum(bs_df$office == office),
                                       scale_sample_size = scale_samp,
                                       data_in = NULL,
                                       verbose = FALSE,
                                       return_program = FALSE,
                                       gc = TRUE)
}
toc(log = TRUE)
rm(blah)
