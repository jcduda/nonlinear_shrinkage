
# STOPPED HERE ADAPT EVERYTHING

# Example fits

library(ggplot2)
library(tidyverse)
library(splines)
library(invgamma)
library(DoseFinding)
library(ToxicR)

load("040_results_simulation_table.RData")
#load("040_results_simulation_raw.RData")
source("./010_functions/020_nlfs_head_fct.R")
source("./010_functions/030_bspline_head_fct.R")
source("./010_functions/010_gen_data_functions.R")
source("./010_functions/011_prep_B.R")
source("./010_functions/021_nlfs_updating_functions.R")
source("./010_functions/040_parametric_bspline_head_function.R")
source("./010_functions/041_parametric_bspline_updating_functions.R")
source("./010_functions/050_p_spline_head_fct.R")

# Case: Truth = hill, n=50 ----------------------------------------------------

res_red <- res_tab %>% filter(n == 50, truth == "hill") %>%
  mutate(method_long = paste0(algorithm,"_",assume,"_",shrinkage)) %>%
  filter(method_long %in% c("parametric_power_NA",
                            "nlfs_power_own_slice",
                            "pspline_NA_NA",
                            "nlfs_hill_own_slice")) %>%
  group_by(method_long) %>%
  filter(row_number() == 1)
#                  algorithm %in% c("parametric_bspline", truth == "hill") %>%
# group_by(assume) %>%
# filter(row_number() == 1)

ids <- as.character(res_red$job.id)

names(results_simulation_raw) <- res_tab$job.id

fits <- results_simulation_raw[ids]
names(fits)

# Redo nlfs(hill) fit to get smooth B-spline grid: #############################
d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "hill")

theta_prior <- get_prior(model = "hill")

set.seed(1234)
res_nlfs_hill <- nlfs(d = d,
                      true_fn = d$true_fn,
                      assume = "hill",
                      shrinkage = "own_slice",
                      n_draw = 12000,
                      n_inner_knots = 15,
                      theta3_mu = theta_prior[1],
                      theta3_sq = theta_prior[2],
                      theta4_mean = theta_prior[3],
                      theta4_var = theta_prior[4],
                      # for assume = hill_power
                      theta5_mu = theta_prior[5],
                      theta5_sq = theta_prior[6])

res_nlfs_hill$B
x_vis <- seq(0, 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_nlfs_hill$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_nlfs_hill$draws$all_beta[, 2001:12000])
nlfs_hill_yhat <- interc_mean + B_vis%*%beta_mean

# Redo P-spline fit to get smooth B-spline grid: ###############################

set.seed(1234)
res_pspline_hill <- pspline(d = d,
                            n_draw = 12000,
                            n_inner_knots = 15)


x_vis <- seq(0, 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_pspline_hill$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_pspline_hill$draws$all_beta[, 2001:12000])
pspline_hill_yhat <- interc_mean + B_vis%*%beta_mean

# Fit parametric Hill-model (oracle scenario) #################################
d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "hill")

theta_prior <- get_prior(model = "hill")

h_log_var <<- log(theta_prior[4] / theta_prior[3]^2 + 1)
h_log_mean <<- log(theta_prior[3]) - h_log_var/2

prior_list <- create_prior_list(normprior(0,1,-100,100), # a or intercept
                                normprior(0,1,-1e4,1e4), # b or Emax
                                normprior(theta_prior[1], sqrt(theta_prior[2]), 0, 100), # c or ED50
                                lnormprior(h_log_mean, sqrt(h_log_var),0,18), # d or h
                                normprior(-1.75, 1)) # for log(sigma^2)

cont_prior <- create_continuous_prior(prior_list, model = "hill", distribution = "normal")

res_param_hill  <- single_continuous_fit(d$x, d$y, prior=cont_prior,
                                         fit_type="mcmc", BMR_TYPE="sd", BMR = 0.1,
                                         burnin = 2000, samples = 10000)

param_hill_y_hat <- rowMeans(predict(res_param_hill, new_doses = x_vis)$Y)

#-------------------------------------------------------------------------------
#
# First Plot: Truth: Hill
#               methods: NLFS(Hill), Param.(Hill), P-Spline
#-------------------------------------------------------------------------------

#pdf("100_fig_example_sim_red_n50_hill_truth.pdf", width = 6, height = 4)
pdf("101_fig_ex_n50_publication.pdf", width = 4, height = 4)
par(mgp = c(2, 1, 0), mar = c(3, 3, 3, 1))
plot(fits[[1]]$x, fits[[1]]$y, pch = 16, col = scales::alpha("black", 0.2),
     xlab = "Dose", ylab = "Posterior mean", main = "Truth = Hill, n = 50")
# pspline
lines(x_vis, pspline_hill_yhat, col = "black", lwd = 2, lty = 1)
#nlfs(Hill)
lines(x_vis, nlfs_hill_yhat, col = "forestgreen", lwd = 3)
#Param(Hill)
lines(x_vis, param_hill_y_hat, col = "magenta", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("NLFS(Hill)", "Param.(Hill)", "P-Spline"),
       col = c("forestgreen", "magenta", "black"),
       lwd = c(2, 3, 2), lty = c(1, 2, 1))
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Redo nlfs(power) fit to get smooth B-spline grid: ############################
d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "hill")

theta_prior <- get_prior(model = "power")

set.seed(1234)
res_nlfs_power <- nlfs(d = d,
                       true_fn = d$true_fn,
                       assume = "power",
                       shrinkage = "own_slice",
                       n_draw = 12000,
                       n_inner_knots = 15,
                       theta3_mu = theta_prior[1],
                       theta3_sq = theta_prior[2],
                       theta4_mean = theta_prior[3],
                       theta4_var = theta_prior[4],
                       # for assume = hill_power
                       theta5_mu = theta_prior[5],
                       theta5_sq = theta_prior[6])

x_vis <- seq(0, 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_nlfs_power$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_nlfs_power$draws$all_beta[, 2001:12000])
nlfs_power_yhat <- interc_mean + B_vis%*%beta_mean


# Redo B-spline fit to get smooth B-spline grid: ###############################

set.seed(1234)
res_bspline_hill <- bspline(d = d,
                            n_draw = 12000,
                            n_inner_knots = 15)


x_vis <- seq(0, 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_bspline_hill$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_bspline_hill$draws$all_beta[, 2001:12000])
bspline_hill_yhat <- interc_mean + B_vis%*%beta_mean

# Parametric(Power) fit ########################################################

d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "power")

theta_prior <- get_prior(model = "power")

h_log_var <<- log(theta_prior[4] / theta_prior[3]^2 + 1)
h_log_mean <<- log(theta_prior[3]) - h_log_var/2

prior_list <- create_prior_list(normprior(0,1,-100,100), # a (intercept)
                                normprior(0,1,-1e4,1e4), # b (scaling)
                                normprior(theta_prior[1], sqrt(theta_prior[2]), 0, 100), # d (power-param)
                                normprior(-1.75, 1)) # for log(sigma^2)
cont_prior <- create_continuous_prior(prior_list, model = "power", distribution = "normal")

res_param_power  <- single_continuous_fit(d$x, d$y, prior=cont_prior,
                                          fit_type="mcmc", BMR_TYPE="sd", BMR = 0.1,
                                          burnin = 2000, samples = 10000)

param_power_yhat <- rowMeans(predict(res_param_power, new_doses = x_vis)$Y)


#-------------------------------------------------------------------------------
#
# Second Plot: Truth: Hill
#               methods: NLFS(Power), B-Spline, Param.(Power), P-Spline
#-------------------------------------------------------------------------------

#pdf("100_fig_example_sim_red_n50_hill_truth.pdf", width = 6, height = 4)
pdf("101_fig_ex_n50_publication_b.pdf", width = 4, height = 4)
par(mgp = c(2, 1, 0), mar = c(3, 3, 3, 1))
plot(fits[[1]]$x, fits[[1]]$y, pch = 16, col = scales::alpha("black", 0.2),
     xlab = "Dose", ylab = "Posterior mean", main = "Truth = Hill, n = 50")
# Param.(power)
lines(x_vis, param_power_yhat, col = "black", lwd = 2, lty = 1)
#nlfs(power)
lines(x_vis, nlfs_power_yhat, col = "forestgreen", lwd = 3)
#B-Spline
lines(x_vis, bspline_hill_yhat, col = "magenta", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("NLFS(power)", "B-Spline",  "Param.(power)"),
       col = c("forestgreen", "magenta", "black"),
       lwd = c(2, 3, 2), lty = c(1, 2, 1))
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# Case: Truth = hilldown, n=50 -------------------------------------------------
#-------------------------------------------------------------------------------

res_red <- res_tab %>% filter(n == 50, truth == "hilldown") %>%
  mutate(method_long = paste0(algorithm,"_",assume,"_",shrinkage)) %>%
  filter(method_long %in% c("parametric_power_NA",
                            "parametric_hill_NA",
                            "nlfs_power_own_slice",
                            "nlfs_hill_own_slice",
                            "nlfs_hill_power_own_slice")) %>%
  group_by(method_long) %>%
  filter(row_number() == 18) %>%
  select(job.id, method_long)


ids <- as.character(res_red$job.id)

names(results_simulation_raw) <- res_tab$job.id

fits <- results_simulation_raw[ids]


# Redo nlfs(Hill+power) fit to get smooth B-spline grid: #######################
d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "hilldown")

theta_prior <- get_prior(model = "hill_power")

set.seed(1234)
res_nlfs_hillpower <- nlfs(d = d,
                           #true_fn = d$true_fn,
                           assume = "hill_power",
                           shrinkage = "own_slice",
                           n_draw = 12000,
                           n_inner_knots = 15,
                           theta3_mu = theta_prior[1],
                           theta3_sq = theta_prior[2],
                           theta4_mean = theta_prior[3],
                           theta4_var = theta_prior[4],
                           # for assume = hill_power
                           theta5_mu = theta_prior[5],
                           theta5_sq = theta_prior[6])


x_vis <- seq(min(d$x), 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_nlfs_hillpower$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_nlfs_hillpower$draws$all_beta[, 2001:12000])
nlfs_hillpower_yhat <- interc_mean + B_vis%*%beta_mean

# Redo P-spline fit to get smooth B-spline grid: ###############################

set.seed(1234)
res_pspline_hilldown <- pspline(d = d,
                                n_draw = 12000,
                                n_inner_knots = 15)


x_vis <- seq(0, 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_pspline_hilldown$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_pspline_hilldown$draws$all_beta[, 2001:12000])
pspline_hilldown_yhat <- interc_mean + B_vis%*%beta_mean


# True mean response ###########################################################

truth_hillpower <- d$true_fn(x_vis)

#-------------------------------------------------------------------------------
#
# Third Plot: Truth: Hill + Downturn
#               methods: NLFS(Hill+Power), P-Spline, Truth
#-------------------------------------------------------------------------------

#pdf("100_fig_example_sim_red_n50_hill_truth.pdf", width = 6, height = 4)
pdf("101_fig_ex_n50_publication_c.pdf", width = 4, height = 4)
par(mgp = c(2, 1, 0), mar = c(3, 3, 3, 1))
plot(fits[[1]]$x, fits[[1]]$y, pch = 16, col = scales::alpha("black", 0.2),
     xlab = "Dose", ylab = "Posterior mean", main = "Truth = Hill + Downturn, n = 50")
# P-Spline
lines(x_vis, pspline_hilldown_yhat, col = "magenta", lwd = 2, lty = 1)
#nlfs(hill+power)
lines(x_vis, nlfs_hillpower_yhat, col = "forestgreen", lwd = 3)
# Truth:
lines(x_vis, truth_hillpower, col = "black", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("NLFS(Hill+power)", "P-Spline",  "Truth"),
       col = c("forestgreen", "magenta", "black"),
       lwd = c(3, 2, 2), lty = c(1, 1, 2))
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



################################################################################

# Forth plot: Compare CIs between param(Hill)+B-spline and nlfs(hill)

################################################################################


# Case: Truth = hill, n=50 ----------------------------------------------------

res_red <- res_tab %>% filter(n == 50, truth == "hill") %>%
  mutate(method_long = paste0(algorithm,"_",assume,"_",shrinkage)) %>%
  filter(method_long %in% c("parametric_power_NA",
                            "nlfs_power_own_slice",
                            "pspline_NA_NA",
                            "nlfs_hill_own_slice")) %>%
  group_by(method_long) %>%
  filter(row_number() == 1)
#                  algorithm %in% c("parametric_bspline", truth == "hill") %>%
# group_by(assume) %>%
# filter(row_number() == 1)

ids <- as.character(res_red$job.id)

names(results_simulation_raw) <- res_tab$job.id

fits <- results_simulation_raw[ids]
names(fits)

# Redo nlfs fit to get smooth B-spline grid: ###################################
d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "hill")

theta_prior <- get_prior(model = "hill")

set.seed(1234)
res_nlfs_hill <- nlfs(d = d,
                      true_fn = d$true_fn,
                      assume = "hill",
                      shrinkage = "own_slice",
                      n_draw = 12000,
                      n_inner_knots = 15,
                      theta3_mu = theta_prior[1],
                      theta3_sq = theta_prior[2],
                      theta4_mean = theta_prior[3],
                      theta4_var = theta_prior[4],
                      # for assume = hill_power
                      theta5_mu = theta_prior[5],
                      theta5_sq = theta_prior[6])


x_vis <- seq(0, 1, 0.01)
B_vis <- prep_B(x_vis, 15, min(d$x), max(d$x), c(0, 1))$B
interc_mean <- mean(res_nlfs_hill$draws$all_interc[2001:12000])
beta_mean <- rowMeans(res_nlfs_hill$draws$all_beta[, 2001:12000])
nlfs_hill_yhat <- interc_mean + B_vis%*%beta_mean

# Get upper and lower limits, but on fine grid:

y_hat_vis <- matrix(res_nlfs_hill$draws$all_interc[2001:12000], byrow = TRUE,
                    nrow = length(x_vis), ncol = length(res_nlfs_hill$draws$all_interc[2001:12000])) +
  B_vis %*% res_nlfs_hill$draws$all_beta[, 2001:12000]

nlfs_upper_vis <- apply(y_hat_vis, 1, function(x) quantile(x, 0.95))
nlfs_lower_vis <- apply(y_hat_vis, 1, function(x) quantile(x, 0.05))

################################################################################

# Param. (Hill) + bspline

d = fits[[1]][c("x", "y")]
d$n <- length(d$x)
d$true_fn <- get_true_fn(truth = "hill")

theta_prior <- get_prior(model = "hill")

set.seed(1234)
res_param_bspline <- parametric_bspline(d = d, n_draw = 12000,
                                        n_inner_knots = 15,
                                        assume = "hill",
                                        theta3_mu = theta_prior[1],
                                        theta3_sq = theta_prior[2],
                                        theta4_mean = theta_prior[3],
                                        theta4_var = theta_prior[4])

# Get the parametric part
samples <- 2001:12000
all_param <- sapply(samples, function(i){
  sigEmax(x_vis, e0 = res_param_bspline$draws$all_interc[i],
          eMax = res_param_bspline$draws$all_emax[i],
          ed50 = res_param_bspline$draws$all_theta[1, i],
          h = res_param_bspline$draws$all_theta[2, i])
})

all_non_param <- B_vis %*% res_param_bspline$draws$all_beta[, samples]

y_hat_vis <- all_param + all_non_param

param_bspline_upper_vis <- apply(y_hat_vis, 1, function(x) quantile(x, 0.95))
param_bspline_lower_vis <- apply(y_hat_vis, 1, function(x) quantile(x, 0.05))

param_bspline_yhat <- rowMeans(y_hat_vis)

# Plot together


pdf("101_fig_ex_n50_publication_d.pdf", width = 4, height = 4)
par(mgp = c(2, 1, 0), mar = c(3, 3, 3, 1))
plot(d$x, d$y, pch = 16, col = scales::alpha("black", 0.3), xlab = "Dose",
     ylab = "Response", main = "Truth = Hill, n = 50")
# nlfs
lines(x_vis, nlfs_hill_yhat, col= "forestgreen", lty = 1, lwd = 2)
polygon(x = c(x_vis, rev(x_vis)), y = c(nlfs_upper_vis, rev(nlfs_lower_vis)),
        col = scales::alpha("forestgreen", 0.2), border = NA)
# param(hill) + B-spline
lines(x_vis, param_bspline_yhat, col= "magenta", lty = 2, lwd = 2)
lines(x_vis, param_bspline_upper_vis, lty = 3, col = "magenta", lwd = 2)
lines(x_vis, param_bspline_lower_vis, lty = 3, col = "magenta", lwd = 2)
legend("bottomright", legend = c("NLFS(Hill): Post. mean ",
                                 "NLFS(Hill): 95% CI",
                                 "Param.(Hill) + B-spline: Post. mean ",
                                 "Param.(Hill) + B-spline: 95% CI"),
       lwd = c(2, 5, 2, 2), col = c("forestgreen", scales::alpha("forestgreen", 0.2),
                                    "magenta", "magenta"),
       lty = c(1, 1, 2, 3), cex = 0.7)
dev.off()
