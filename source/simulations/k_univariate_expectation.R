####################################################################
# Julia Wrobel
# August 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# comparing fperm to CSR to Kinhom,
####################################################################

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
source(here::here("source", "simulate_scSpatialSim.R"))
source(here::here("source", "utils_k.R"))

###############################################################
## set simulation design elements
###############################################################

n = c(1000, 2000, 5000, 10000)
abundance = c(0.01, 0.1, 0.2)
type = c("hom", "inhom", "homClust", "inhomClust")
nperm = 1000
seed_start = 1000
N_iter = 50

params = expand.grid(seed_start = seed_start,
                     type = type,
                     n = n,
                     abundance = abundance)  %>%
  mutate(m = n * abundance) %>%
  filter(m >=5)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("output", "univariate_expectation", "varyAbundance_scSim"), Date), showWarnings = FALSE)


## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 4
  #scenario = 3
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
abundance = params$abundance[scenario]
type = params$type[scenario]
SEED.START = params$seed_start[scenario]

results = vector("list", length = N_iter)
for(iter in 1:N_iter){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter
  set.seed(seed.iter)

  # simulate data
  if(type %in% c("hom", "inhom")){
    ppp_obj <- mxsim_univariate(n, abundance, type)
  }else{
    ppp_obj <- sim_scSpatial(n, abundance, type)
  }


  ################################################################################
  ##
  # Calculate Ripley's K and fperm statistics
  k_full = get_k(ppp_obj, nperm = nperm)

  lambda_n = n
  lambda_m = m
  res = mutate(k_full, n = ppp_obj$n, m = subset(ppp_obj, marks == "immune")$n) %>%
    mutate(iter = iter,
           scenario = scenario,
           seed = seed.iter,
           type = type,
           nperm = nperm,
           lambda_n = lambda_n,
           abundance = abundance)


  results[[iter]] = res

} # end for loop


filename = paste0(here::here("output", "univariate_expectation", "varyAbundance_scSim", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


