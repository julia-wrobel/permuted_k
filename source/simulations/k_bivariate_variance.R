####################################################################
# Julia Wrobel
# October 2024
#
# This file produces simulations for bivariate K under different data generation mechanisms
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
suppressPackageStartupMessages(library(scSpatialSIM))



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
source(here::here("source", "utils_k_bivariate.R"))
source(here::here("source", "simulate_scSpatialSim.R"))
source(here::here("source", "get_permutation_distribution.R"))

###############################################################
## set simulation design elements
###############################################################

n = c(1000, 2000, 5000, 10000)
abundance = c(0.01, 0.1, 0.2)
type = c("hom", "inhom", "homClust", "inhomClust")
nperm = 1000
seed_start = 1000
N_iter = 1000
maxiter = (seq(1, N_iter, by = 100)-1) + 100

params = expand.grid(seed_start = seed_start,
                     type = type,
                     n = n,
                     abundance = abundance,
                     maxiter = (seq(1, N_iter, by = 100)-1) + 100) %>%
  mutate(m = n * abundance) %>%
  filter(m >=5)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("output", "bivariate_variance", "varyAbundance"), Date), showWarnings = FALSE)

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
  if(type %in% c("hom", "inhom")){
    ppp_obj <- mxsim(n, abundance, type, bivariate = TRUE)
  }else{
    ppp_obj <- sim_scSpatial(n, abundance, type, bivariate = TRUE)
  }

  ################################################################################
  ##
  # Calculate Ripley's K and fperm statistics
  k_kamp = get_k_power_biv(ppp_obj)
  k_perm = get_k_power_permOnly_biv(ppp_obj, nperm = nperm)

  lambda_n = n
  lambda_m = m
  res = mutate(bind_rows(k_kamp, k_perm),n = ppp_obj$n,
               m1 = subset(ppp_obj, marks == "immune1")$n,
               m2 = subset(ppp_obj, marks == "immune2")$n
               ) %>%
    mutate(iter = iter_vec[i],
           scenario = scenario,
           seed = seed.iter,
           type = type,
           lambda_n = lambda_n,
           abundance = abundance)


  results[[i]] = res
} # end for loop


filename = paste0(here::here("output", "bivariate_variance", "varyAbundance", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


