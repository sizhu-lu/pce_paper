# Replication files for ``Principal Stratification with Continuous Post-Treatment Variables''
# PART II: Generate Monte Carlo samples under three regimes

# load the required packages
library(np)
library(ks)
library(pracma)
library(parallel)
library(abind)
library(doParallel)
library(dplyr)
library(plyr)

# load all functions from part I
source("part_I_all_functions.R")

# ------------------------------------------------------------------------ #
# Regime 1: both treatment probability and outcome mean models are correct #
# ------------------------------------------------------------------------ #

# estimation function for one sample
# rho_both_correct
get_results_one_mc_sample <- function(n=500, n_divisions=100, suffix=0, rho=0) {
  # generate data
  n = n
  p = 2
  X = data.frame(matrix(rnorm(n*p, 0, 1), n, p))
  # tilde_X = ((X + 0.25) ^ 2 - 1) / sqrt(2)
  
  # tp model
  Z = rbinom(n, 1, 1/2)
  # tp_mis <- apply(tilde_X, 1, sum) / 2
  # Z = rbinom(n, 1, plogis(tp_mis))
  # ps model
  e_rho <- rnorm(n)
  rho <- rho
  S1 = 0.6 * X[,1] + 0.8 * (sqrt(rho) * e_rho + sqrt(1-rho) * rnorm(n))
  S0 = 0.6 * X[,2] + 0.8 * (sqrt(rho) * e_rho + sqrt(1-rho) * rnorm(n))
  
  # outcome model
  Y1 = S1 - X[,1] + rnorm(n)
  Y0 = S0 - X[,2] + rnorm(n)
  # Y1 = S1 + 2 * sqrt(2) * tilde_X[, 1] + rnorm(n)
  # Y0 = S0 + 2 * sqrt(2) * tilde_X[, 2] + rnorm(n)
  
  # observed data
  S = Z * S1 + (1-Z) * S0
  Y = Z * Y1 + (1-Z) * Y0
  data = list(X=X, Z=Z, S=S, Y=Y)
  df = data.frame(data)
  
  if (nrow(X) == 0) {stop("X data is empty!")}
  
  # estimate
  result <- point_estimator(Z, X, S, Y, n_divisions, copula_type='gaussian', rho=rho)
  save(result, file = paste0("simulation_results/rho_all_correct/rho_", rho, "_all_correct_n_", n, "_mc_", suffix, "_20231227.RData"))
}

# get the number of available cores and set
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
# stopCluster(cl)

# simulation by parallel computing
MC <- 500
message("Start time: ", Sys.time())
for (rho in c(0, 0.2, 0.5)) {
  for (n in c(200, 500, 1000, 2000)) {
    foreach(i=1:MC) %dopar% {
      get_results_one_mc_sample(n=n, n_divisions=100, suffix=i, rho=rho)
    }
  }
}
beepr::beep()

# --------------------------------------------------------- #
# Regime 2: the treatment probability model is misspecified #
# --------------------------------------------------------- #

# estimation function for one sample
# rho_tp_wrong
get_results_one_mc_sample <- function(n=500, n_divisions=100, suffix="est", rho=0) {
  
  # generate data
  n = n
  p = 2
  X = data.frame(matrix(rnorm(n*p, 0, 1), n, p))
  tilde_X = ((X + 0.25) ^ 2 - 1) / sqrt(2)
  
  # tp model
  # Z = rbinom(n, 1, 1/2)
  tp_mis <- apply(tilde_X, 1, sum) / 2
  Z = rbinom(n, 1, plogis(tp_mis))
  # ps model
  e_rho <- rnorm(n)
  rho <- rho
  S1 = 0.6 * X[,1] + 0.8 * (sqrt(rho) * e_rho + sqrt(1-rho) * rnorm(n))
  S0 = 0.6 * X[,2] + 0.8 * (sqrt(rho) * e_rho + sqrt(1-rho) * rnorm(n))
  
  # outcome model
  Y1 = S1 - X[,1] + rnorm(n)
  Y0 = S0 - X[,2] + rnorm(n)
  # Y1 = S1 + 2 * sqrt(2) * tilde_X[, 1] + rnorm(n)
  # Y0 = S0 + 2 * sqrt(2) * tilde_X[, 2] + rnorm(n)
  
  # observed data
  S = Z * S1 + (1-Z) * S0
  Y = Z * Y1 + (1-Z) * Y0
  data = list(X=X, Z=Z, S=S, Y=Y)
  df = data.frame(data)
  
  if (nrow(X) == 0) {stop("X data is empty!")}
  
  # estimate
  result <- point_estimator(Z, X, S, Y, n_divisions, copula_type='gaussian', rho=rho)
  save(result, file = paste0("simulation_results/rho_tp_wrong/rho_", rho, "_tp_wrong_n_", n, "_mc_", suffix, "_20231227.RData"))
}

# get the number of available cores and set
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

# simulation by parallel computing
MC <- 500
message("Start time: ", Sys.time())
for (rho in c(0, 0.2, 0.5)) {
  for (n in c(200, 500, 1000, 2000)) {
    foreach(i=1:MC) %dopar% {
      get_results_one_mc_sample(n=n, n_divisions=100, suffix=i, rho=rho)
    }
  }
}
beepr::beep()

# -------------------------------------------------- #
# Regime 3: the outcome mean models are misspecified #
# -------------------------------------------------- #

# estimation function for one sample
# rho_om_wrong
get_results_one_mc_sample <- function(n=500, n_divisions=100, suffix="est", rho=0) {
  
  # generate data
  n = n
  p = 2
  X = data.frame(matrix(rnorm(n*p, 0, 1), n, p))
  tilde_X = ((X + 0.25) ^ 2 - 1) / sqrt(2)
  
  # tp model
  Z = rbinom(n, 1, 1/2)
  # tp_mis <- apply(tilde_X, 1, sum) / 2
  # Z = rbinom(n, 1, plogis(tp_mis))
  # ps model
  e_rho <- rnorm(n)
  rho <- rho
  S1 = 0.6 * X[,1] + 0.8 * (sqrt(rho) * e_rho + sqrt(1-rho) * rnorm(n))
  S0 = 0.6 * X[,2] + 0.8 * (sqrt(rho) * e_rho + sqrt(1-rho) * rnorm(n))
  
  # outcome model
  # Y1 = S1 + X[,1] + rnorm(n)
  # Y0 = S0 + X[,2] + rnorm(n)
  Y1 = S1 - 2 * sqrt(2) * tilde_X[, 1] + rnorm(n)
  Y0 = S0 - 2 * sqrt(2) * tilde_X[, 2] + rnorm(n)
  
  # observed data
  S = Z * S1 + (1-Z) * S0
  Y = Z * Y1 + (1-Z) * Y0
  data = list(X=X, Z=Z, S=S, Y=Y)
  df = data.frame(data)
  
  if (nrow(X) == 0) {stop("X data is empty!")}
  
  # estimate
  result <- point_estimator(Z, X, S, Y, n_divisions, copula_type='gaussian', rho=rho)
  save(result, file = paste0("simulation_results/rho_om_wrong/rho_", rho, "_om_wrong_n_", n, "_mc_", suffix, "_20231227.RData"))
}

# get the number of available cores and set
no_cores <- detectCores() - 2
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

# simulation by parallel computing
MC <- 500
message("Start time: ", Sys.time())
for (rho in c(0, 0.2, 0.5)) {
  for (n in c(200, 500, 1000, 2000)) {
    foreach(i=1:MC) %dopar% {
      get_results_one_mc_sample(n=n, n_divisions=100, suffix=i, rho=rho)
    }
  }
}
beepr::beep()
