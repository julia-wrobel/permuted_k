

# function to combine files that are being stored as a list, store this in utils.
# reads in a file name, the file should be an Rda object that is a list of data frames that have the same structure
# and can be just bound together
merge_files = function(file, variance = FALSE){
  load(file)
  res = map_dfr(results, bind_rows)

  if(variance){
    res = res %>%
      filter(r == 0.15, holes == FALSE)

  }else{
    res = res %>%
      mutate(type_char = case_when(
        type == 1 ~ "hom",
        type == 2 ~ "inhom",
        type == 3 ~ "homClust",
        TRUE ~ "inhomClust"
      ))
  }


  res

}




