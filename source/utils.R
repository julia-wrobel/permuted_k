

# function to combine files that are being stored as a list, store this in utils.
# reads in a file name, the file should be an Rda object that is a list of data frames that have the same structure
# and can be just bound together
merge_files = function(file, variance = FALSE){
  load(file)
  res = map_dfr(results, bind_rows)

  if(variance){
    res = res #%>%
      #filter(r == 0.15, holes == FALSE)
  }

  res
}



# merge and process files from real data analysis
merge_files_vectraPolaris = function(file, variance = FALSE, bivariate = FALSE){
  if(bivariate){
    load(here::here("output", "vpData_bivariate", file))
  }else{
    load(here::here("output", "vpData", file))
  }

  sample_id = str_replace(file, "RDA", "im3")

  if(variance){
    res = results$variance

  }else{
    res = results$expectation
  }

  res = res %>% mutate(sample_id = sample_id)
  res
}



# function to extract data from survival analysis
get_coxph = function(dat, rval, doc_method = "kamp", get_resi = TRUE){
  dat = filter(dat, r == rval, method == doc_method)  %>%
    mutate(doc = doc / 1000)

  mod = coxph(Surv(survival_time, event) ~ age + stage + p_immune + doc,
              data = dat)

  if(get_resi){
    resi_object = resi(mod, data = dat)

    tidy(mod, exp = TRUE) %>%
      mutate(r = rval, method = doc_method,
             resi = coef(resi_object)$RESI,
             lower = coef(resi_object)$"2.5%",
             upper = coef(resi_object)$"97.5%")
  }else{
    tidy(mod, exp = TRUE) %>%
      mutate(r = rval, method = doc_method)
  }
}



