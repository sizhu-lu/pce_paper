# Replication files for ``Principal Stratification with Continuous Post-Treatment Variables''
# PART IV: Application codes and generate Figure 1

# load required packages
library(dplyr)
library(doParallel)
library(ggplot2)
library(tidyr)
library(latex2exp)
library(labeling)

# load all functions from the file
source("part_I_all_functions.R")

# load the simulated data
load("application_simulated_data.RData")

# post-treatment variable: d_inc1 (change in expected monthly income 6 months after the experiment)
# outcome variable: d_cc0 (change in average monthly spending after the survey)
df <- simulated_data
Z <- df$treat
X <- df %>% select(female, age, safety_score, w)
S <- df$d_inc1
Y <- df$d_cc0

# get point estimators
rho_list <- c(0, 0.1, 0.2, 0.5, 0.8)
for (rho_ind in 1:length(rho_list)) {
  rho = rho_list[rho_ind]
  point_est <- point_estimator(Z, X, S, Y, n_divisions=100, copula_type='gaussian', rho=rho)
  save(point_est, file = paste0("application_results/point_res_rho_", rho, ".RData"))
}

## nonparametric bootstrap for variance estimation
# get the number of available cores and set
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores, type="FORK")
registerDoParallel(cl)

n_boot = 200
n_divisions = 100
copula_type='gaussian'
message("Start time: ", Sys.time())
for (rho_ind in 1:length(rho_list)) {
  rho = rho_list[rho_ind]
  foreach (i = 1:n_boot) %dopar% {
    set.seed(i+2023)
    # nonparametric bootstrap
    n <- length(Z)
    X <- as.matrix(X)
    id_boot = sample(1:n, n, replace = TRUE)
    est_boot = point_estimator(Z[id_boot], X[id_boot, ], S[id_boot], Y[id_boot], n_divisions, copula_type, rho,
                               weighting_function_vectorized, g_function_vectorized)
    save(est_boot, file = paste0("application_boot_res/boot_rho_", rho, "_est_", i, ".RData"))
  }
}
message("End time: ", Sys.time())
beepr::beep()

# summarize the results
for (rho in rho_list) {
  all_boot_res <- c()
  for (i in 1:n_boot) {
    load(paste0("application_boot_res/boot_rho_", rho, "_est_", i, ".RData"))
    all_boot_res = rbind(all_boot_res, c(est_boot, rho))
    colnames(all_boot_res) <- c("eif_beta1", "eif_beta0", "eif_alpha", 
                                "tp_pd_beta1", "tp_pd_beta0", "tp_pd_alpha", 
                                "pd_om_beta1", "pd_om_beta0", "pd_om_alpha", "rho")
    write.csv(all_boot_res, file = paste0("application_results/application_boot_all_res_rho", rho, ".csv"))
  }
}

# replication codes for the barplot in Figure 1
rho_list <- c(0, 0.1, 0.2, 0.5, 0.8)
res <- c()
for (rho in rho_list){
  load(paste0("application_results/point_res_rho_", rho, ".RData"))
  all_boot_res <- read.csv(paste0("application_results/application_boot_all_res_rho", rho, ".csv"))
  boot_ci <- c()
  for (i in 2:10) {
    boot_ci <- rbind(boot_ci, c(quantile(all_boot_res[, i], 0.025), quantile(all_boot_res[, i], 0.975)))
  }
  res <- rbind(res, cbind(data.frame(cbind(point_est, boot_ci)), 
                          rho = paste0("rho:", rho),
                          est = c(rep("eif", 3), rep("tp+pd", 3), rep("pd+om", 3)),
                          coef = rep(c("beta[1]", "beta[0]", "alpha"), 3)))
}
colnames(res) <- c("point", "lb", "ub", "rho", "est", "coef")
df <- data.frame(res)

# barplot
bar_plot <- ggplot(df, aes(x = factor(est, level=c("eif", "tp+pd", "pd+om")), y = point)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymax = ub, ymin = lb),
                width=0.2, size=0.2) +
  facet_grid(factor(coef, levels=c("beta[1]", "beta[0]", 'alpha')) ~ rho,
             scales = "free_y",
             labeller = label_parsed) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank()) +
  xlab("estimator")
print(bar_plot)

