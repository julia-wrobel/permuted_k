####################################################################
# Julia Wrobel
# August 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# focusing on the variance/power. Runs 1000 iterations in chunks of 50 at a time.
####################################################################

#suppressPackageStartupMessages()

suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(spatstat.random))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(spatstat.explore))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(scSpatialSIM))
suppressPackageStartupMessages(library(survival))



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
source(here::here("source", "simulate_scSpatialSim.R"))


###############################################################
## set simulation design elements
###############################################################

n = c(5000)
abundance = c(0.1)
type = c("inhomClust")
beta_val = c(0, 0.1, 0.5, 2)
rho = c(-0.5, .5) # correlation of covariates
seed_start = 2000
N_iter = 1000
n_subj = c(100, 500, 1000)
maxiter = (seq(1, N_iter, by = 100)-1) + 100

params = expand.grid(seed_start = seed_start,
                     type = type,
                     n = n,
                     abundance = abundance,
                     beta_val = beta_val,
                     rho = rho,
                     n_subj = n_subj,
                     maxiter = (seq(1, N_iter, by = 100)-1) + 100) %>%
  mutate(m = n * abundance) %>%
  filter(m >=5)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("output", "univariate_survival"), Date), showWarnings = FALSE)

## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 52
  it = 100
  n_subj = 500
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
  it = 100
}


###############################################################
## start simulation code
###############################################################

###############################################################
## set simulation design elements
###############################################################
n = params$n[scenario]
abundance = params$abundance[scenario]
n_subj = params$n_subj[scenario]
m = n * abundance
type = params$type[scenario]
beta_val = params$beta_val[scenario]
rho = params$rho[scenario]
SEED.START = params$seed_start[scenario]
maxiter = params$maxiter[scenario]

iter_vec = (maxiter-99):maxiter

results = vector("list", length = it)
for(i in 1:it){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter_vec[i]
  set.seed(seed.iter)


  ################################################################################
  ##
  # Calculate Ripley's K and fperm statistics
  kamp = rnorm(n_subj)
  k = kamp + rnorm(n_subj, mean = 7)
  r = 1
  id = 1:n_subj


  kvals = data.frame(id = id, r = r, k = k, kamp = kamp)


  # simulate survival data and fit cox models
  fits = get_coxPH(n_subj, beta_val, rho, kvals)


  fits = fits %>%
    mutate(iter = iter_vec[i],
           scenario = scenario,
           seed = seed.iter,
           type = type,
           abundance = abundance,
           rho = rho,
           n_subj = n_subj,
           beta_val = beta_val)


  results[[i]] = fits
} # end for loop


filename = paste0(here::here("output", "univariate_survival", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


