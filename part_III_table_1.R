# Replication files for ``Principal Stratification with Continuous Post-Treatment Variables''
# PART III: Process results and generate Table 1

# load the required packages
library(dplyr)
library(plyr)
library(reshape)
library(pracma)

# function to compute the ground truth
ground_truth <- function(rho) {
  eta1 <- inv(matrix(c(1, rho * 0.8 ^ 2, 0, rho * 0.8 ^ 2, 1, 0, 0, 0, 1), nrow = 3)) %*% c(1.6, rho * 0.8 ^ 2, 0)
  eta1_mis <- inv(matrix(c(1, rho * 0.8 ^ 2, 0, rho * 0.8 ^ 2, 1, 0, 0, 0, 1), nrow = 3)) %*% c(1.6, rho * 0.8 ^ 2, 0.125)
  return(eta1[1] - eta1[2])
}

ground_truth_res <- c()
for (rho in c(0, 0.2, 0.5)) {
  ground_truth_res <- rbind(ground_truth_res, c(rho, ground_truth(rho)))
}

# summarize all results
MC <- 500
result_table <- c()
for (rho in c(0, 0.2, 0.5)) {
  sub_table <- c()
  # all_correct
  for (n in c(200, 500, 1000, 2000)) {
    all_result <- c()
    for (i in 1:MC) {
      file <- paste0("simulation_results/rho_all_correct/rho_", rho, "_all_correct_n_", n, "_mc_", i, "_20231227.RData")
      load(file)
      centered_result <- result - c(ground_truth(rho), -ground_truth(rho), 0)
      all_result <- rbind(all_result, centered_result)
    }
    sub_table <- rbind(sub_table, c(rho, n, apply(all_result, 2, mean), 
                                    apply(all_result, 2, sd), 
                                    sqrt(apply(all_result^2, 2, mean))))
  }
  # tp_wrong
  for (n in c(200, 500, 1000, 2000)) {
    all_result <- c()
    for (i in 1:MC) {
      file <- paste0("simulation_results/rho_tp_wrong/rho_", rho, "_tp_wrong_n_", n, "_mc_", i, "_20231227.RData")
      load(file)
      centered_result <- result - c(ground_truth(rho), -ground_truth(rho), 0)
      all_result <- rbind(all_result, centered_result)
    }
    sub_table <- rbind(sub_table, c(rho, n, apply(all_result, 2, mean), 
                                    apply(all_result, 2, sd), 
                                    sqrt(apply(all_result^2, 2, mean))))
  }
  # om_wrong
  for (n in c(200, 500, 1000, 2000)) {
    all_result <- c()
    for (i in 1:MC) {
      file <- paste0("simulation_results/rho_om_wrong/rho_", rho, "_om_wrong_n_", n, "_mc_", i, "_20231227.RData")
      load(file)
      centered_result <- result - c(ground_truth(rho), -ground_truth(rho), 0.125)
      all_result <- rbind(all_result, centered_result)
    }
    sub_table <- rbind(sub_table, c(rho, n, apply(all_result, 2, mean), 
                                    apply(all_result, 2, sd), 
                                    sqrt(apply(all_result^2, 2, mean))))
  }
  
  # store the results
  result_table <- rbind(result_table, sub_table)
}

# generate Table 1
ind_list <- seq(1:9) * 3
beta_1_table <- result_table[, c(1, 2, ind_list)][, c(1, 2, 3, 6, 9, 4, 7, 10, 5, 8, 11)]
beta_0_table <- result_table[, c(1, 2, ind_list+1)][, c(1, 2, 3, 6, 9, 4, 7, 10, 5, 8, 11)]
alpha_table <- result_table[, c(1, 2, ind_list+2)][, c(1, 2, 3, 6, 9, 4, 7, 10, 5, 8, 11)]
combined_table <- data.frame(rbind(beta_1_table, beta_0_table, alpha_table))
colnames(combined_table) <- c("rho", "n", 
                              "eif_bias", "eif_sd", "eif_rmse", 
                              "tp+pd_bias", "tp+pd_sd", "tp+pd_rmse",
                              "pd+om_bias", "pd+om_sd", "pd+om_rmse")
combined_table$regime <- rep(c(rep("all_correct", 4), rep("tp_wrong", 4), rep("om_wrong", 4)), 9)
combined_table$estimator <- c(rep("beta1", 36), rep("beta0", 36), rep("alpha", 36))
# the first 36 rows in combined_table are estimators on beta1 and are collected in Table 1
