
################################################################################
# Code for Figure 5: Testosterone example
################################################################################


library(dplyr)
library(DoseFinding)
library(splines)
library(invgamma)

source("../02_simulation_study/010_functions/020_nlfs_head_fct.R")
source("../02_simulation_study/010_functions/030_bspline_head_fct.R")
source("../02_simulation_study/010_functions/010_gen_data_functions.R")
source("../02_simulation_study/010_functions/011_prep_B.R")
source("../02_simulation_study/010_functions/021_nlfs_updating_functions.R")
source("../02_simulation_study/010_functions/040_parametric_bspline_head_function.R")
source("../02_simulation_study/010_functions/041_parametric_bspline_updating_functions.R")
source("../02_simulation_study/010_functions/050_p_spline_head_fct.R")

d1 <- readxl::read_xlsx("testosterone_data.xlsx")
d2 <- d1 %>% dplyr::select("log10(TT Observed + 1)", "Age Observed")
colnames(d2) <- c("log10TTp1", "age")
d2$TT <- 10^(d2$log10TTp1) - 1
d <- d2 %>% # filter(age <= 60) %>%
  filter(age <= 85) %>%
  mutate(x = age, y = TT) %>%
  dplyr::select(x, y)

d <- as.list(d)
d$n <- length(d$x)


# Load NLFS(Hill) fit result

load("./020_testosterone_nlfs_res_50knots.RData")

# The literature model: --------------------------------------------------------

a <- 0.04655
b <- -0.05311
c0 <- 0.05123
d0 <- -0.00793
e <- -0.01222
f0 <- 0.00058
g <- 0.00069

their_fn <- function(x) 10^((a + c0*x + e* x^2 + g*x^3) /
                              (1 + b * x +  d0*x^2 + f0*x^3)) - 1

# P-Spline ---------------------------------------------------------------------

# Takes a minute
set.seed(1234)
res_pspline <- pspline(d = d,
                            n_draw = 2000,
                            n_inner_knots = 50,
                       boundary_knots = c(min(d$x), max(d$x)))


################################################################################
# Figure 5: NLFS, Literature model, P-spline
################################################################################

pdf("fig_050_testosterone.pdf", width = 5, height = 4)
par(mgp = c(2, 1, 0), mar = c(3.5, 3, 2, 1))
# mean Posterior
plot(d$x, d$y, col = scales::alpha("black", 0.05), xlab = "Age [years]",
     ylab = "Total Testosterone [nmol/L]")
# P-Spline
lines(d$x, rowMeans(res_pspline$draws$all_y_hat[,1001:2000]), col = "orange", lwd = 2)
# NLFS
lines(d$x, rowMeans(res_nlfs$draws$all_y_hat[,1001:2000]), col = "forestgreen", lwd = 2)
# Their Curve
curve(their_fn, col = "magenta", add = T, lwd = 2, lty = 2)

legend("topright", legend = c("NFLS(Hill)", "Literature model",  "P-spline"),
       col = c("forestgreen", "magenta", "orange"), lwd = 2, lty = c(1, 2, 1))

dev.off()

