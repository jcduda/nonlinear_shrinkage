

# Code for Figure 1 and 2 in the manuscript

# Figure 1 #####################################################################

library(invgamma)

# We mathematically derived what is (prop. to) the density of omega:

pdf("fig_01_omega_prior_labels_with_Omega0.pdf", width = 6, height = 4)
par(mar = c(2.5, 3, 2, 0.6), mgp = c(1.5, 0.5, 0), cex = 1.2)
a = b = 0.5
f_omega <- function(x) 1/2 * x^(a-1) * (1-x)^(b-1)
curve(f_omega, from = 0.01, to = 0.99,
      ylab = "Density up to constants", xaxt = "n",
      xlab = expression(omega==1/(1+tau^2)))
axis(side=1, at=c(0, 1))
text(0.8, 4.5, expression(Shrinkage))
text(0.8, 4, expression((f~'in'~{Omega[0]^{Theta}})))

text(0.2, 4.5, expression(No~shrinkage))
text(0.2, 4, expression((f~outside~Omega[0]^{Theta})))
dev.off()

# Figure 2 #####################################################################


source("../02_simulation_study/010_functions/010_gen_data_functions.R")

# stopped here: fix: use actual E0 intercept!

all_non_lin_theta <- c(ED50 = 0.3, h = 6)
hill1 <- function(x) sigEmax(x, e0 = 0.2, eMax = 1,
                             ed50 = all_non_lin_theta["ED50"],
                             h = all_non_lin_theta["h"])



pdf("fig_02_hill_schematic_curve.pdf", width = 4, height = 3)
par(mar = c(3, 3, 0.3, 0.1), mgp = c(2, 1, 0))
curve(hill1,  xlab = "Dose", ylab = "Response", ylim = c(0, 1.2),
      yaxt = "n")

axis(2, at = c(0, 0.2, 0.7, 1.2), labels = c(0, 0.2, 0.7, 1.2))

# Theta 1
arrows(x0 = 0, y0 = 0, y1= 0.2, length = 0.1)
arrows(x0 = 0, y0 = 0.2, y1= 0, length = 0.1)
text(x=0.05, y = 0.1, label = expression(theta[1]))

# Theta 2
abline(h=1.2, lty = 2, col = scales::alpha("black", 0.5))
arrows(x0 = 0, y0 = 0.2, y1= 1.2, length = 0.1)
arrows(x0 = 0, y0 = 1.2, y1= 0.2, length = 0.1)
text(x=0.03, y = 0.9, label = expression(theta[2]))

# Theta 3
lines(x = c(0, 0.3), y = c(0.7, 0.7), lty = 2, col = scales::alpha("black", 0.5))
lines(x = c(0.3, 0.3), y = c(0.7, 0.1), lty = 2, col = scales::alpha("black", 0.5))
text(x=0.3, y = 0.05, label = expression(theta[3]))

# Theta 4
arrows(x0= 0.35, x1= 0.42, y0 = 0.5, y1=0.8, length = 0.1)
text(x=0.41, y = 0.61, label = expression(theta[4]))
dev.off()

