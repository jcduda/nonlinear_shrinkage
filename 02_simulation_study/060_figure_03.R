
################################################################################
# Code for Figure 3 (a & b): Boxplots of summarized, reduced simulation results
################################################################################

library(ggplot2)
library(tidyverse)
library(patchwork)
load("040_results_simulation_table.RData")

#############
# Figure 3a
#############

# Only the case where the truth is the Hill model
curr_truth <- "hill"

res_reduced <- res_tab %>%
  filter(truth == curr_truth,
         (shrinkage %in% c("own_slice", NA) ) | (algorithm == "parametric_bspline" & shrinkage == "half_cauchy"),
         assume %in% c("hill", "power", "hill_power", NA)) %>%
  # To take a look at the relevant combination/filtering settings
  # dplyr::select(true_fn, sigma_sq, method, assume, prior, shrinkage) %>% distinct() %>%
  mutate(method_long = paste0(algorithm,"_", assume)) %>%
  mutate(method_long = gsub("_NA", "", method_long)) %>%
  mutate(method_long = factor(method_long,
                              levels = c("parametric_power", "nlfs_power", "bspline",
                                         "parametric_bspline_power", "pspline",
                                         "parametric_bspline_hill", "nlfs_hill_power", "nlfs_hill", "parametric_hill")),
         n = paste0("n = ", n),
         n = factor(n, levels = paste0("n = ", c(50, 100, 200, 500)))) %>%
  # Filter further for a reduced graphics
  filter(method_long %in% c("parametric_power", "nlfs_power", "pspline", "nlfs_hill", "parametric_hill"),
         n == "n = 50")

# Help function to allow line break in labels
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

p1 <- res_reduced %>%
  ggplot(., aes(fill = method_long, y = RMSE, x = method_long)) +
  geom_boxplot(fill="white") +
  theme_bw() +
  scale_fill_manual(values="white") +
  scale_x_discrete(labels = addline_format(c("Param. (misspec.)", "NLFS (misspec.)",
                                             "P-Spline", "NLFS (correct)", "Param. (correct)")))+
  labs(fill = "Method", x = "Method", title = "a) Truth: Hill (n=50)")  +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position = "none",# "bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1), colour = "black"),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.title = element_text(size = rel(1.2)))

p1


#############
# Figure 3b
#############

# Only the case where the truth is a Hill model with a downturn effect
curr_truth <- "hilldown"

res_reduced <- res_tab %>%
  filter(truth == curr_truth,
         (shrinkage %in% c("own_slice", NA) ) | (algorithm == "parametric_bspline" & shrinkage == "half_cauchy"),
         assume %in% c("hill", "power", "hill_power", NA)) %>%

  mutate(method_long = paste0(algorithm,"_", assume)) %>%
  mutate(method_long = gsub("_NA", "", method_long)) %>%
  mutate(method_long = factor(method_long,
                              levels = c("nlfs_power", "parametric_bspline_power", "nlfs_hill",
                                         "parametric_bspline_hill", "nlfs_hill_power",
                                         "parametric_power",  "bspline",
                                         "pspline", "parametric_hill")),
         n = paste0("n = ", n),
         n = factor(n, levels = paste0("n = ", c(50, 100, 200, 500)))) %>%
  # Filter further to reduce graphic
  filter(method_long %in% c("nlfs_power", "parametric_bspline_power", "nlfs_hill",
                            "parametric_bspline_hill", "nlfs_hill_power"),
         n == "n = 50")


p2 <- res_reduced %>%
  ggplot(., aes(fill = method_long, y = RMSE, x = method_long)) +
  geom_boxplot(fill="white") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c("NLFS (power)", "Param.(power) +B-spline ",
                                             "NLFS (Hill)", "Param.(Hill) +B-spline", "NLFS (Hill+power)"))) +
  labs(fill = "Method", x = "Method", title = " b) Truth: Hill + Downturn (n=50)")  +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position = "none",
        legend.title=element_blank(),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1), colour = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = rel(1.2)))


# Put 3a (p1) and 3b (p2) next to each other, share legend

p1 + p2 + plot_layout(guides = "collect", widths = c(1, 1.2))

ggsave("fig_03_sim_res_reduced.pdf", width = 10, height = 5)
