####################################################################
# Julia Wrobel
# August 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# focusing on the variance/power. Runs 1000 iterations in chunks of 50 at a time.
####################################################################

#suppressPackageStartupMessages()

suppressPackageStartupMessages(library(spatstat.random))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tictoc))


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

n = c(100, 500, 2000, 5000)
abundance = c(0.01, 0.1, 0.5)
type = c("hom", "inhom", "homClust",  "inhomClust", "inhomTightClust")

seed_start = 1000
N_iter = 1000
nperm = c(100, 10000)
maxiter = (seq(1, N_iter, by = 100)-1) + 100

params = expand.grid(seed_start = seed_start,
                     type = type,
                     nperm = nperm,
                     n = n,
                     abundance = abundance,
                     maxiter = (seq(1, N_iter, by = 100)-1) + 100) %>%
  mutate(m = n * abundance) %>%
  filter(m >5)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("output", "univariate_variance", "varyAbundance_permOnly"), Date), showWarnings = FALSE)

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
abundance = params$abundance[scenario]
m = n * abundance
nperm = params$nperm[scenario]

type = params$type[scenario]
SEED.START = params$seed_start[scenario]
maxiter = params$maxiter[scenario]

iter_vec = (maxiter-99):maxiter

results = vector("list", length = 100)
for(i in 1:100){
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
  k_full = get_k_power_permOnly(ppp_obj$full, nperm = nperm)


  lambda_n = n
  lambda_m = m
  res = mutate(k_full, holes = FALSE, n = ppp_obj$full$n, m = subset(ppp_obj$full, marks == "immune")$n) %>%
    mutate(iter = iter_vec[i],
           scenario = scenario,
           seed = seed.iter,
           type = type,
           lambda_n = lambda_n,
           abundance = abundance,
           nperm = nperm)


  results[[i]] = res
} # end for loop


filename = paste0(here::here("output", "univariate_variance", "varyAbundance_permOnly", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


