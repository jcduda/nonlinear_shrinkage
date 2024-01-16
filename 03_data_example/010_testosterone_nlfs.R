# Model testosterone data set as application example: Now: below60

# Set working directory to 03_data_example
# setwd("./03_data_example")

library(dplyr)
library(DoseFinding)
library(splines)
library(invgamma)

d1 <- readxl::read_xlsx("testosterone_data.xlsx")

source("../02_simulation_study/010_functions/010_gen_data_functions.R")
source("../02_simulation_study/010_functions/011_prep_B.R")
source("../02_simulation_study/010_functions/020_nlfs_head_fct.R")
source("../02_simulation_study/010_functions/021_nlfs_updating_functions.R")

# Prepare data set -------------------------------------------------------------

d2 <- d1 %>% dplyr::select("log10(TT Observed + 1)", "Age Observed")

colnames(d2) <- c("log10TTp1", "age")

d2$TT <- 10^(d2$log10TTp1) - 1


d <- d2 %>% filter(age <= 85) %>%
  mutate(x = age, y = TT) %>%
  dplyr::select(x, y)

d <- as.list(d)
d$n <- length(d$x)

# Running NLFS can take long and requires ~ 50 GB memory

theta_start <- c(15, 10)


set.seed(1234)
res_nlfs = nlfs(n_draw = 2000, n_inner_knots = 50,
                d = d,
                assume = "hill",
                sigma_sq_half_cauchy = TRUE,
                theta3_mu = 15,
                theta3_sq = 4,
                theta4_mean = 10,
                theta4_var = 5,
                theta_start = c(15, 10),
                boundary_knots = c(min(d$x), max(d$x)))

save(res_nlfs, file = "020_testosterone_nlfs_res_50knots.RData")
