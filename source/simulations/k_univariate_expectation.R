####################################################################
# Julia Wrobel
# August 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# comparing fperm to CSR to Kinhom,
####################################################################

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
source(here::here("source", "utils_fast.R"))

###############################################################
## set simulation design elements
###############################################################

n = c(500, 2000) # we will need to go bigger
nm = c(50, 100, 500)
type = c("hom", "inhom", "homClust", "inhomClust")

seed_start = 1000
N_iter = 50

params = expand.grid(seed_start = seed_start,
                     type = type,
                     n = n,
                     nm = nm)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path("./output", Date), showWarnings = FALSE)

## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 1
  N_iter = 2
}else{
  # defined from batch script params
  scenario <- commandArgs(trailingOnly=TRUE)
}

###############################################################
## start simulation code
###############################################################


###############################################################
## set simulation design elements
###############################################################
n = params$n[scenario]
nm = params$nm[scenario]
type = params$type[scenario]
SEED.START = params$seed_start[scenario]

results = vector("list", length = N_iter)
for(iter in 1:N_iter){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter
  set.seed(seed.iter)

  # simulate data
  ppp_obj <- NULL
  attempt <- 1
  while(is.null(ppp_obj) && attempt <= 5) {
    attempt <- attempt + 1
    try(
      ppp_obj <- mxsim_univariate(n, nm, type)
    )
  }

  par = c(iter, scenario, seed.iter, type, ppp_obj$full$n, subset(ppp_obj$full, marks == "immune")$n,
          ppp_obj$holes$n, subset(ppp_obj$holes, marks == "immune")$n, n, nm)

  par = matrix(par, nrow = 1)
  ################################################################################
  ##
  # Calculate Ripley's K and fperm statistics
  k_full = get_k(ppp_obj$full)
  k_holes = get_k(ppp_obj$holes)


  results_mat = cbind(k_full, k_holes)

  n1 = c("iter","scenario", "seed", "type", "n", "nm", "n_hole", "nm_hole", "lambda_n", "lambda_nm")
  n2 = colnames(k_full)
  n3 = paste0(colnames(k_holes), "_hole")

  colnames(results_mat) <- c(n2, n3)
  colnames(par) <- n1

  results[[iter]] = bind_cols(as_tibble(results_mat), as_tibble(par))
} # end for loop


filename = paste0(here::here("output", Date), "/k_univariate_", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


