

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
merge_files_vectraPolaris = function(file, variance = FALSE){
  load(here::here("output", "vpData", file))
  sample_id = str_replace(file, "RDA", "im3")

  if(variance){
    res = results$variance

  }else{
    res = results$expectation
  }

  res = res %>% mutate(sample_id = sample_id)
  res
}




