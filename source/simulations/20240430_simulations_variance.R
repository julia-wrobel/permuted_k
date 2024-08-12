####################################################################
# Julia Wrobel
# April 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# and calculates null variance for p-value
####################################################################

library(spatstat.random)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(tictoc)


# get parallelization set up
library(foreach)
library(doParallel)
n_cores = detectCores() - 4
registerDoParallel(n_cores)


### Notes

* https://arxiv.org/pdf/1006.1567
* https://www.researchgate.net/publication/236331582_A_Statistical_Test_for_Ripley's_K_Function_Rejection_of_Poisson_Null_Hypothesis

* Permutation test can be approcimated by a normal distribution with the same mean and variance
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4682210/

###############################################################
## define or source functions used in code below
###############################################################
source(here::here("source", "simulate_ppp.R"))
source(here::here("source", "utils_fast.R"))

###############################################################
## set simulation design elements
###############################################################

#n = c(500, 2000) # we will need to go bigger
n = 2000
nm = c(20, 50)
type = c("hom", "inhom", "homClust", "inhomClust")

seed_start = 1000
N_iter = 500


params = expand.grid(seed_start = seed_start,
										 type = type,
										 n = n,
										 nm = nm)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())


###############################################################
## start simulation code
###############################################################


## need to parallelize this part
## actually, parallelize the iterations?
for(scenario in 1:nrow(params)){

	###############################################################
	## set simulation design elements
	###############################################################
	n = params$n[scenario]
	nm = params$nm[scenario]
	type = params$type[scenario]
	SEED.START = params$seed_start[scenario]



	k_univariate = foreach(iter = 1:N_iter, .combine = 'rbind') %dopar% {

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

	filename = paste0(here::here("output", "simulation_results", Date), "/k_univariate", scenario, ".RDA")
	save(k_univariate,
			 file = filename)
} # end scenario loop

###############################################################
## end sim
###############################################################


