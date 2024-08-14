library(tidyverse)
library(stringr)
library(viridis)
library(spatstat.random)
library(patchwork)

files = list.files(here::here("output", "20240813"),
                   full.names =  TRUE)




# function to combine files that are being stored as a list, store this in utils.
# reads in a file name, the file should be an Rda object that is a list of data frames that have the same structure
# and can be just bound together
merge_files = function(file){
  load(file)
  map_dfr(results, bind_rows) %>%
    # clean parameter scenarios
    mutate(type_char = case_when(
      type == 1 ~ "hom",
      type == 2 ~ "inhom",
      type == 3 ~ "homClust",
      TRUE ~ "inhomClust"
    ))
}

expectation_df = map_dfr(files, merge_files)


# rearrange dataset
expectation_df %>%
  select(contains("time_"), starts_with("k"), everything()) %>%
  pivot_longer(time_kinhom:kfperm_thin_hole,
               names_to = "method", values_to = "value") %>%
  mutate(measurement = ifelse(grepl("time_", method), "time", "value"),
         method = str_remove(method, "time_")) %>%
  pivot_wider(names_from = measurement, values_from = value) %>%
  relocate(method, value, time)
  mutate(holes = ifelse(grepl("_hole", method), TRUE, FALSE),
         method = str_remove(method, "_hole"))


df = expectation_df %>%
  mutate(k_theo = khat - ktheo,
         k_inhom = kinhomhat - kinhomtheo,
         k_perm = khat - kperm,
         k_fperm = khat - kfperm,
         k_fpermThin = khat - kfperm_thin) %>%
  select(r, scenario, type_char, n, nm, lambda_n, lambda_nm,
         contains("k_"), contains("time"))
,,x
