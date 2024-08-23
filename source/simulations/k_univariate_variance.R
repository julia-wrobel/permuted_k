####################################################################
# Julia Wrobel
# August 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# focusing on the variance/power. Runs 1000 iterations in chunks of 50 at a time.
####################################################################

#suppressPackageStartupMessages()

library(spatstat.random)
library(spatstat.geom)
library(spatstat.explore)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(tictoc)


wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}


###############################################################
## define or source functions used in code below
###############################################################
source(here::here("source", "simulate_ppp.R"))
source(here::here("source", "utils_k.R"))
source(here::here("source", "get_permutation_distribution.R"))

###############################################################
## set simulation design elements
###############################################################

n = c(500, 2000)
m = c(50, 100, 500)
type = c("hom", "inhom", "homClust", "inhomClust")

seed_start = 1000
N_iter = 1000
maxiter = (seq(1, N_iter, by = 50)-1) + 50

params = expand.grid(seed_start = seed_start,
                     type = type,
                     n = n,
                     m = m,
                     maxiter = (seq(1, N_iter, by = 50)-1) + 50)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("output", "univariate_variance"), Date), showWarnings = FALSE)

## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 1
  N_iter = 2
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}


###############################################################
## start simulation code
###############################################################

###############################################################
## set simulation design elements
###############################################################
n = params$n[scenario]
m = params$m[scenario]
type = params$type[scenario]
SEED.START = params$seed_start[scenario]
maxiter = params$maxiter[scenario]

iter_vec = (maxiter-49):maxiter

results = vector("list", length = 50)
for(i in 1:50){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter_vec[i]
  set.seed(seed.iter)

  # simulate data
  ppp_obj <- NULL
  attempt <- 1
  while(is.null(ppp_obj) && attempt <= 5) {
    attempt <- attempt + 1
    try(
      ppp_obj <- mxsim_univariate(n, m, type)
    )
  }

  ################################################################################
  ##
  # Calculate Ripley's K and fperm statistics
  k_full = get_k_power(ppp_obj$full)
  k_holes = get_k_power(ppp_obj$holes)


  res = bind_rows(mutate(k_full, holes = FALSE, n = ppp_obj$full$n, m = subset(ppp_obj$full, marks == "immune")$n),
            mutate(k_holes, holes = TRUE, n = ppp_obj$holes$n, m = subset(ppp_obj$holes, marks == "immune")$n)) %>%
    mutate(iter = iter_vec[i],
           scenario = scenario,
           seed = seed.iter,
           type = type,
           lambda_n = n,
           lambda_m = m)


  results[[i]] = res
} # end for loop


filename = paste0(here::here("output", "univariate_variance", Date), "/",scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


