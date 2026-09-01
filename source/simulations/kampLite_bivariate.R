####################################################################
# Julia Wrobel
#
# This file produces simulations for KAMP lite under different thinning proportions
# and data generation scenarios, for the bivariate (immune1 x immune2) setting.
# Varies over different abundances for the two bivariate cell types.
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
source(here::here("source", "utils_k_bivariate.R"))
source(here::here("source", "get_permutation_distribution.R"))

###############################################################
## set simulation design elements
###############################################################

n = c(5000, 10000, 20000)
abundance1 = c(0.001, 0.05, 0.01)
abundance2 = c(0.01, 0.1, 0.2)
type = c("inhom", "inhomClust")
seed_start = 1000
N_iter = 50

params = expand.grid(seed_start = seed_start,
                     type = type,
                     abundance1 = abundance1,
                     abundance2 = abundance2,
                     n = n)  %>%
  mutate(m1 = n * abundance1,
         m2 = n * abundance2) %>%
  filter(m1 >=3, m2 >=3)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("output", "kamplite_bivariate"), Date), showWarnings = FALSE)


## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 3
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
m1 = params$m1[scenario]
m2 = params$m2[scenario]
abundance1 = params$abundance1[scenario]
abundance2 = params$abundance2[scenario]
type = params$type[scenario]
thinning_proportions = c(0, 0.25, 0.5, 0.75, 0.9)
SEED.START = params$seed_start[scenario]

results = vector("list", length = N_iter)
for(iter in 1:N_iter){
  # set seed
  seed.iter = (SEED.START - 1)*N_iter + iter
  set.seed(seed.iter)

  # simulate data
  abundance = c(abundance1, abundance2)
  if(type %in% c("hom", "inhom")){
    ppp_obj <- mxsim(n, abundance, type, bivariate = TRUE)
  }else{
    ppp_obj <- sim_scSpatial(n, abundance, type, bivariate = TRUE)
  }


  ################################################################################
  ##
  ## go over all the thinning options you are interested in
  k_lite = map_dfr(thinning_proportions, get_kamplite_biv,
          ppp_obj = ppp_obj, rvec = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2))


  lambda_n = n
  lambda_m1 = m1
  lambda_m2 = m2
  res = mutate(k_lite, n = ppp_obj$n,
              m1 = subset(ppp_obj, marks == "immune1")$n,
              m2 = subset(ppp_obj, marks == "immune2")$n) %>%
    mutate(iter = iter,
           scenario = scenario,
           seed = seed.iter,
           type = type,
           lambda_n = lambda_n,
           abundance1 = abundance1,
           abundance2 = abundance2)


  results[[iter]] = res

} # end for loop


filename = paste0(here::here("output", "kamplite_bivariate", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################

