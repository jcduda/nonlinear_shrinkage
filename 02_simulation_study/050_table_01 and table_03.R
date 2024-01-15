
###########################################################################
# Generate LaTeX code for Table 1 & 3 in manuscript
###########################################################################


library(tidyverse)
library(xtable)
load("040_results_simulation_table.RData")

# Table 1: Method overview ##############################################################
res_tab %>% select(algorithm, assume, shrinkage) %>% distinct %>%
  mutate(shrinkage = ifelse(shrinkage == "half_cauchy", "a=b=0.5", "a=0.5, b=exp(-klog(n)/2)")) %>%
  xtable()


# Table 3: Simulation results ###########################################################

res_tab %>%
  mutate(datagen_long = paste0("n_", n, "_",truth),
         method_long = paste0(algorithm,"_",assume,"_",shrinkage)) %>%
  mutate(method_long = gsub("half_cauchy", "HC", method_long),
         method_long = gsub("own_slice", "OS", method_long),
         method_long = gsub("_NA", "", method_long),
         method_long = gsub("parametric", "param", method_long)) %>%
  dplyr::select(datagen_long, method_long, RMSE) %>%
  group_by(datagen_long, method_long) %>%
  summarise(RMSE_mean = mean(RMSE, na.rm = T),
            RMSE_SD = sqrt(var(RMSE, na.rm=T))) %>%
  mutate(RMSE_res = paste0(round(RMSE_mean, 3), " (", round(RMSE_SD, 3), ")")) %>%
  ungroup() %>%
  select(-c("RMSE_mean", "RMSE_SD")) %>%
  mutate(datagen_long = factor(datagen_long, levels = c("n_50_hill",
                                                        "n_50_power",
                                                        "n_50_hilldown",
                                                        "n_100_hill",
                                                        "n_100_power",
                                                        "n_100_hilldown",
                                                        "n_200_hill",
                                                        "n_200_power",
                                                        "n_200_hilldown",
                                                        "n_500_hill",
                                                        "n_500_power",
                                                        "n_500_hilldown"))) %>%
  arrange(by = datagen_long) %>%

  pivot_wider(id_cols = method_long, names_from = datagen_long, values_from = RMSE_res) %>%

  mutate(method_long = factor(method_long, levels =
                                c("nlfs_hill_OS", "nlfs_power_OS", "nlfs_hill_power_OS",
                                  "nlfs_hill_HC", "nlfs_power_HC", "nlfs_hill_power_HC",
                                  "param_hill", "param_power",
                                  "bspline", "pspline",
                                  "param_bspline_hill_HC", "param_bspline_power_HC"
                              ))) %>%
  arrange(method_long) %>%

  xtable()
