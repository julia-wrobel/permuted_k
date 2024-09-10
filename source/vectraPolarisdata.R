####################################################################
# Julia Wrobel
# August 2024
#
# This file produces simulations for univariate K under different data generation mechanisms
# focusing on the variance/power. Runs 1000 iterations in chunks of 50 at a time.
####################################################################

#suppressPackageStartupMessages()


suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(tidyverse))
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
#source(here::here("source", "utils_k.R"))
#source(here::here("source", "get_permutation_distribution.R"))

###############################################################
## set simulation design elements
###############################################################

params = expand.grid(seed_start = seed_start,
                     type = type,
                     n = n,
                     m = m,
                     maxiter = (seq(1, N_iter, by = 500)-1) + 500)

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())


## define number of simulations and parameter scenario
if(doLocal) {

  sample_index = 1

  # load data from VectraPolarisData package and process here
  suppressPackageStartupMessages(library(VectraPolarisData))
  spe_ovarian <- HumanOvarianCancerVP()

  ## Assays slots
  assays_slot <- assays(spe_ovarian)
  intensities_df <- assays_slot$intensities
  nucleus_intensities_df<- assays_slot$nucleus_intensities
  rownames(nucleus_intensities_df) <- paste0("nucleus_", rownames(nucleus_intensities_df))
  membrane_intensities_df<- assays_slot$membrane_intensities
  rownames(membrane_intensities_df) <- paste0("membrane_", rownames(membrane_intensities_df))

  # colData and spatialData
  colData_df <- colData(spe_ovarian)
  spatialCoords_df <- spatialCoords(spe_ovarian)

  # clinical data
  patient_level_ovarian <- metadata(spe_ovarian)$clinical_data %>%
    # create binary stage variable
    dplyr::mutate(stage_bin = ifelse(stage %in% c("1", "2"), 0, 1))

  ovarian <- as.data.frame(cbind(colData_df,
                            spatialCoords_df,
                            t(intensities_df),
                            t(nucleus_intensities_df),
                            t(membrane_intensities_df))
  ) %>%
    dplyr::rename(cd19 = cd19_opal_480,
                  cd68 = cd68_opal_520,
                  cd3 = cd3_opal_540,
                  cd8 = cd8_opal_650,
                  ier3 = ier3_opal_620,
                  pstat3 = p_stat3_opal_570,
                  ck = ck_opal_780,
                  ki67 = ki67_opal_690,
                  x = cell_x_position,
                  y = cell_y_position) %>%
    dplyr::select(contains("id"), tissue_category, x,y, contains("phenotype"), ck:dapi) %>%
    # remove control subjects
    dplyr::filter(sample_id %in% patient_level_ovarian$sample_id) %>%
    mutate(immune = ifelse(phenotype_cd19 == "CD19+" | phenotype_cd8 == "CD8+" |
                             phenotype_cd3 == "CD3+" | phenotype_cd68 == "CD68+", "immune", "background"),
           immune = factor(immune, levels = c("immune", "background"))) %>%
    filter(tissue_category == "Tumor") %>%
    select(cell_id, sample_id, x, y, immune)

  rm(spe_ovarian, assays_slot, intensities_df, nucleus_intensities_df, membrane_intensities_df, colData_df, spatialCoords_df)


  # save processed data
  save(ovarian, file = here::here("data", "processed_ovarian_data.Rda"))

}else{
  sample_index <- as.numeric(commandArgs(trailingOnly=TRUE))
  # load preprocessed data

  load(file = here::here("data", "processed_ovarian_data.Rda"))
}


###############################################################
## set which sample to run
###############################################################


# subset data to a single sample
ids = unique(ovarian$sample_id)
ovarian = ovarian %>%
  filter(sample_id == ids[sample_index])


marksvar = "immune"
rvalues = seq(0, 200, length.out = 200) # fine grid of R values

w = convexhull.xy(ovarian[["x"]], ovarian[["y"]])

# define pp_obj
pp_obj = ppp(ovarian[["x"]], ovarian[["y"]], window = w, marks = ovarian[[marksvar]])


################################################################################
## Expectation
################################################################################

k_expectation = get_k(pp_obj, rvec = rvalues)


k_expectation = as_tibble(k_expectation)

################################################################################
## Variance
################################################################################

k_variance = get_k_power(pp_obj, rvec = rvalues)

# variance code not currently optimized for multiple radii- would be good to implement this
k_variance = k_variance %>%
  mutate(time = ifelse(method == "kepd", time/length(rvalues), time))

################################################################################
## Save results
################################################################################

results = list(expectation = k_expectation,
               variance = k_variance)

filename = paste0(here::here("output", "VPdata"), "/", str_remove(ids[sample_index], ".im3"), ".RDA")
save(results,
     file = filename)

###############################################################
## end analysis
###############################################################


