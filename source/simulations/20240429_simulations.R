####################################################################
# Julia Wrobel
# April 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# comparing fperm to CSR to Kinhom,
####################################################################

library(spatstat.random)
library(spatstat.explore)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(tictoc)


# get parallelization set up
library(foreach)
library(doParallel)
n_cores = detectCores() - 2
registerDoParallel(n_cores)

###############################################################
## define or source functions used in code below
###############################################################
source(here::here("source", "simulate_ppp.R"))
source(here::here("source", "utils_fast.R"))

###############################################################
## set simulation design elements
###############################################################

#n = c(500, 2000) # we will need to go bigger
# expected number of background cells
n = 2000

# expected number of immune cells
nm = c(50)
# clustering scenario for the simulated data. 'hom' is homogeneous background, homogeneous immune cells.
  # 'inhom' is inhomogeneous background, 'homClust' is homogeneous background, clustered immune cells,
  # 'inhomClust' is inhomogeneous background, clustered immune cells
type = c("hom", "inhom", "homClust", "inhomClust")

seed_start = 1000

# number of simulated datasets generated. Chose small number just to test the code
N_iter = 10


params = expand.grid(seed_start = seed_start,
										 type = type,
										 n = n,
										 nm = nm)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())


###############################################################
## start simulation code
###############################################################


for(scenario in 1:nrow(params)){

	###############################################################
	## set simulation design elements
	###############################################################
	n = params$n[scenario]
	nm = params$nm[scenario]
	type = params$type[scenario]
	SEED.START = params$seed_start[scenario]


 # parallelize each simulated dataset run
	k_univariate = foreach(iter = 1:N_iter, .combine = 'rbind') %dopar% {

	  # set seed
	  seed.iter = (SEED.START - 1)*N_iter + iter
	  set.seed(seed.iter)

	  # simulate data
	  # sometimes the cluster simulation scenario breaks, which is why I'm doing error catching here
	  ppp_obj <- NULL
	  attempt <- 1
	  while(is.null(ppp_obj) && attempt <= 5) {
	    attempt <- attempt + 1
	    try(
	      ppp_obj <- mxsim_univariate(n, nm, type)
	    )
	  }


	  par = c(iter, scenario, seed.iter, type, ppp_obj$full$n, subset(ppp_obj$full, marks == "immune")$n, n, nm)

	  par = matrix(par, nrow = 1)
		################################################################################
		##
    # Calculate Ripley's K and fperm statistics
	  k_full = get_k(ppp_obj$full)




	  results_mat = cbind(matrix(par, nrow = 3, ncol = length(par), byrow = T), k_full)

	  n1 = c("iter","scenario", "seed", "type", "n", "nm", "lambda_n", "lambda_nm")
	  n2 = colnames(k_full)

		colnames(results_mat) <- c(n1, n2)

	  results_mat
	} # end N_iter foreach loop

	filename = paste0(here::here("output", "simulation_results", Date), "_k_univariate", scenario, ".RDA")
	save(k_univariate,
			 file = filename)
} # end scenario loop

###############################################################
## end sim
###############################################################


