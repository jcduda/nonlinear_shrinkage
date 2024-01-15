library(batchtools)
library(data.table)


# Run the simulation using the batchtools R-package for job scheduling on an HPC
# See: https://github.com/mllg/batchtools
# You might need to setup a configuration file appropriate for your HPC
# You can also run the simulation locally, without the parallel HPC scheduling,
# but it will take very long.

################################################################################

# Creaty a registry. Load all required source files and packages.

reg = makeExperimentRegistry(file.dir = "./registry",
                             packages = c("data.table", "dplyr", "tidyr", "DoseFinding",
                                          "invgamma", "splines", "ToxicR"),
                             source = c("./010_functions/010_gen_data_functions.R",
                                       "./010_functions/011_prep_B.R",
                                       "./010_functions/012_reduce_res.R",
                                       "./010_functions/020_nlfs_head_fct.R",
                                       "./010_functions/021_nlfs_updating_functions.R",
                                       "./010_functions/030_bspline_head_fct.R",
                                       "./010_functions/040_parametric_bspline_head_function.R",
                                       "./010_functions/041_parametric_bspline_updating_functions.R",
                                       "./010_functions/050_p_spline_head_fct.R"),
                             seed = 1)

##########################################
# Stochastic data generation as problem  #
##########################################


generate_data <- function(data, job, truth, sigma_sq, n, ...){

  stoch_data <- gen_data(true_fn = get_true_fn(truth),
                         sigma_sq = sigma_sq,
                         n = n)

  list(stoch_data = stoch_data)
}

addProblem(name = "generate_data", data = NA, fun = generate_data, seed = 42)



################################################################################
# Algorithms                                                                   #
################################################################################

# The following algorithms are considered (cf. Table 1 in manuscript)

# 1. nlfs
# 2. parametric
# 3. bspline
# 4. parametric_bspline
# 5. pspline


# 1. nlfs ----------------------------------------------------------------------

nlfs.wrapper = function(data, job, instance,
                        assume, shrinkage, n_draw, n_inner_knots, burnin){

  theta_prior <- get_prior(model = assume)
  d = instance$stoch_data

  res <- tryCatch({nlfs(d = d,
                true_fn = d$true_fn,
                assume = assume,
                shrinkage = shrinkage,
                n_draw = n_draw,
                n_inner_knots = n_inner_knots,
                theta3_mu = theta_prior[1],
                theta3_sq = theta_prior[2],
                theta4_mean = theta_prior[3],
                theta4_var = theta_prior[4],
                # for assume = hill_power
                # (For the case of shrinkage into combined function spaces)
                theta5_mu = theta_prior[5],
                theta5_sq = theta_prior[6])},
                error = function(cond) return(NA))

  res_red <- reduce_res(res = res, method = "nlfs", true_fn = d$true_fn,
                        samples = (burnin+1):(n_draw))
  return(res_red)
}

# 2. parametric-----------------------------------------------------------------

parametric.wrapper = function(data, job, instance, assume, n_draw, burnin){

  d = instance$stoch_data
  theta_prior <<- get_prior(model = assume)
  print(theta_prior)

  if(assume == "hill"){
    h_log_var <<- log(theta_prior[4] / theta_prior[3]^2 + 1)
    h_log_mean <<- log(theta_prior[3]) - h_log_var/2

    prior_list <- create_prior_list(normprior(0,1,-100,100), # a or intercept
                               normprior(0,1,-1e4,1e4), # b or Emax
                               normprior(theta_prior[1], sqrt(theta_prior[2]), 0, 100), # c or ED50
                               lnormprior(h_log_mean, sqrt(h_log_var),0,18), # d or h
                               normprior(-1.75, 1)) # for log(sigma^2)
  }

  if(assume == "power"){
    prior_list <- create_prior_list(normprior(0,1,-100,100), # a (intercept)
                               normprior(0,1,-1e4,1e4), # b (scaling)
                               normprior(theta_prior[1], sqrt(theta_prior[2]), 0, 100), # d (power-param)
                               normprior(-1.75, 1)) # for log(sigma^2)
  }

  cont_prior <- create_continuous_prior(prior_list, model = assume, distribution = "normal")

  res  <- single_continuous_fit(d$x, d$y, prior=cont_prior,
                                fit_type="mcmc", BMR_TYPE="sd", BMR = 0.1,
                                burnin = burnin, samples = n_draw - burnin)

  res_red <- reduce_res(res = res, method = "parametric", true_fn = d$true_fn)

  return(res_red)
}

# 3. bspline -------------------------------------------------------------------

bspline.wrapper = function(data, job, instance,
                       n_draw, n_inner_knots, burnin){


  d = instance$stoch_data

  res <- tryCatch({bspline(d=d, n_draw = n_draw,
                           n_inner_knots = n_inner_knots)},

                  error = function(cond) return(NA))

  res_red <- reduce_res(res = res, method = "bspline", true_fn = d$true_fn,
                        samples = (burnin+1):(n_draw))
  return(res_red)
}

# 4. parametric_bspline --------------------------------------------------------

parametric_bspline.wrapper = function(data, job, instance, assume, shrinkage,
                                      n_draw, n_inner_knots, burnin){

  theta_prior <- get_prior(model = assume)
  d = instance$stoch_data

  res <- tryCatch({parametric_bspline(d=d, n_draw = n_draw,
                                      n_inner_knots = n_inner_knots,
                                      assume = assume,
                                      shrinkage = shrinkage,
                                      true_fn = d$true_fn,
                                      theta3_mu = theta_prior[1],
                                      theta3_sq = theta_prior[2],
                                      theta4_mean = theta_prior[3],
                                      theta4_var = theta_prior[4])},
                  error = function(cond) return(NA))

  res_red <- reduce_res(res = res, method = "parametric_bspline", true_fn = d$true_fn,
                        samples = (burnin+1):(n_draw))
  return(res_red)
}

# 5. pspline -------------------------------------------------------------------

pspline.wrapper = function(data, job, instance,
                           n_draw, n_inner_knots, burnin){
  d = instance$stoch_data

  res <- tryCatch({pspline(d=d, n_draw = n_draw,
                           n_inner_knots = n_inner_knots)},

                  error = function(cond) return(NA))

  res_red <- reduce_res(res = res, method = "pspline", true_fn = d$true_fn,
                        samples = (burnin+1):(n_draw))
  return(res_red)
}

################################################################################
# Problem design
################################################################################

# data generation (problem design for batchtools) is factorial
pdes = list(generate_data = CJ(truth = c("hill", "hilldown", "power"),
                               sigma_sq = 0.005,
                               n = c(50, 100, 200, 500)
)
)

################################################################################
# Add all algorithms
################################################################################

# 1. nlfs
# 2. parametric
# 3. bspline
# 4. parametric_bspline
# 5. pspline


addAlgorithm(name = "nlfs", fun = nlfs.wrapper)
addAlgorithm(name = "parametric", fun = parametric.wrapper)
addAlgorithm(name = "bspline", fun = bspline.wrapper)
addAlgorithm(name = "parametric_bspline", fun = parametric_bspline.wrapper)
addAlgorithm(name = "pspline", fun = pspline.wrapper)

####################
# algorithm design #
####################

all_assume = c("hill", "power")
all_shrinkage = c("own_slice", "half_cauchy")
n_draw = 12000
burnin = 2000
n_inner_knots = 15



ades = list(
  # Only nlfs can also assume hill_power (i.e.: shrinkage into combined function spaces)
  nlfs = CJ(assume = c(all_assume, "hill_power"), shrinkage = all_shrinkage, n_draw = n_draw, burnin = burnin, n_inner_knots = n_inner_knots),
  parametric = CJ(assume = all_assume, n_draw = n_draw, burnin = burnin),
  bspline = CJ(n_draw = n_draw, burnin = burnin, n_inner_knots = n_inner_knots),
  pspline = CJ(n_draw = n_draw, burnin = burnin, n_inner_knots = n_inner_knots),
  parametric_bspline = CJ(assume = all_assume, shrinkage = "half_cauchy", n_draw = n_draw, burnin = burnin, n_inner_knots = n_inner_knots)
)


###################
# Add experiments #
###################

# 1000 replicates per simulation scenario

addExperiments(pdes, ades, repls = 1000)

# Take a look:
summarizeExperiments(by = c("problem", "algorithm", "assume", "shrinkage"))

###########################################################
# Optional: testing single jobs before submitting all jobs
###########################################################


# easy scenario: Clear time-dependency, little noise, large n
id = head(findExperiments(algo.name = "parametric_bspline",
                          prob.pars = (truth == "hill" &
                                       n == 50)),
          1)

testJob(id= id)

# find added experiments:

added_jobs <- findExperiments(algo.name = "nlfs",
                            algo.pars = (assume == "hill_power"))
all_jobs <- findJobs()

################################################################################
# Group jobs into chunks to not send too many single jobs
################################################################################


all_jobs[,chunk:=chunk(all_jobs$job.id, n.chunks = 1000)] # ~145 jobs per chunk



################################################################################
# Submit all jobs
################################################################################

# Memory and walltime requirements might differ on different HPCs
submitJobs(all_jobs, resources = list(walltime = 28800L, memory = 8192))


# check if jobs are done
getStatus()

################################################################################
# If all jobs are done:
################################################################################

# Summarize results into one object (might take a while)
results_simulation_raw <- reduceResultsList()
# Saved bject not shared in Git as it is too big. Can be provided upon request.
save(results_simulation_raw, file = "030_results_simulation_raw.RData")
# extract RMSE as performance measure
RMSEs <- lapply(results_simulation_raw, function(x) x$RMSE) %>% unlist()

# Prepare to add simulation settings to each performance result
job_RMSE_done <- data.frame(job.id = findDone(), RMSE = RMSEs)
pars = unwrap(getJobPars())
res_tab = ijoin(pars, job_RMSE_done)
save(res_tab, file = "040_results_simulation_table")


