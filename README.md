# Replication package for Principal Stratification with Continuous Post-Treatment Variables: Nonparametric Identification and Semiparametric Estimation (Lu Jiang Ding 2023)

The replication folder contains the following files and directories: 
  four R scripts for simulation and application studies, 
  two directories storing our raw results for simulation and application studies, 
  an RData file containing simulated data for the application study. The real data is not publicly available so we generate simulated data to illustrate our methods. 

## 1. Replication of Table 1

The replication process for Table 1 consists of three R scripts:

(1) "part_I_all_functions.R" contains all necessary functions for the simulation study.

(2) "part_II_generate_mc_results.R" contains our data-generating processes for all three regimes and produces raw simulation results. 

(3) "part_III_table_1.R" summarizes raw results into a table containing all information required in Table 1. 
We include our raw results from part II in the directories "rho_all_correct/", "rho_tp_wrong/", and "rho_om_wrong/". Readers may also run part III using our raw results.

## 2. Replication of Figure 1

Figure 1 in the paper is based on real data from _Yin, X. (2022). Learning in the limit: Income inference from credit extension_. SSRN4254400. 
The real data is not publicly available, therefore, we generate simulated data from the original data. The script "part_IV_figure_1.R" can be run using the simulated data; however, the resulting barplot will differ from Figure 1 in the paper due to the differences between the simulated and real data. To replicate Figure 1 exactly, readers may use the raw results stored in "application_results/" and execute the final section of "part_IV_figure_1.R".
