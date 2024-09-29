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

###############################################################
## set simulation design elements
###############################################################

n = c(100, 500, 2000, 5000)
abundance = c(0.01, 0.1, 0.5)
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
dir.create(file.path(here::here("output", "univariate_expectation", "varyAbundance"), Date), showWarnings = FALSE)


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
abundance = params$abundance[scenario]
m = n * abundance
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
      ppp_obj <- mxsim_univariate(n, m, type)
    )
  }

  par = c(iter, scenario, seed.iter, type, ppp_obj$full$n, subset(ppp_obj$full, marks == "immune")$n, n, m, abundance)

  par = matrix(par, nrow = 1)
  ################################################################################
  ##
  # Calculate Ripley's K and fperm statistics
  k_full = get_k(ppp_obj$full, nperm = nperm)



  results_mat = k_full

  n1 = c("iter","scenario", "seed", "type", "n", "m", "lambda_n", "lambda_m", "abundance")
  n2 = colnames(k_full)

  colnames(results_mat) <- c(n2)
  colnames(par) <- n1

  results[[iter]] = bind_cols(as_tibble(results_mat), as_tibble(par))
} # end for loop


filename = paste0(here::here("output", "univariate_expectation", "varyAbundance", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


