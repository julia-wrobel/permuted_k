

# define function for simulating univariate mIF data. Returns object with and without holes
## lambda_n: intensity for background cells
## lambda_m: intensity for marker positive cells
## holes: should an image be simulated with or without holes
## type: defines the distribution of the point process- homogeneous, inhomogeneous, or clustered
sim_scSpatial <- function(lambda_n,
                          abundance,
                          type = c("homClust", "inhomClust"),
                          bivariate = FALSE){


  wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

  if(bivariate){
    sim_object = CreateSimulationObject(sims = 1, cell_types = 2, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    sim_object = GenerateCellPositivity(sim_object,
                                        k = 25,
                                        sdmin = .7, sdmax = .71,
                                        step_size = 0.1, cores = 1,
                                        probs = c(0.0, 1, 1))

    if(type == "inhomClust"){
      sim_object = GenerateHoles(sim_object, step_size = 0.1, cores = 1)
    }

    # get dataframe with the info you want
    pp = CreateSpatialList(sim_object, single_df = TRUE) %>%
      rename(immune1 = `Cell 1 Assignment`,
             immune2 = `Cell 2 Assignment`)

    phat1 = sum(pp$immune1)/nrow(pp)
    phat2 = sum(pp$immune2)/nrow(pp)

    if(phat1 > abundance){
      nhat = nrow(pp)
      nthin = round((phat1 - abundance) * nhat)
      indices = which(pp$immune1 == 1)
      pp$immune1[sample(indices, nthin)] <- 0
    }
    if(phat2 > abundance){
      nhat = nrow(pp)
      nthin = round((phat2 - abundance) * nhat)
      indices = which(pp$immune2 == 1)
      pp$immune2[sample(indices, nthin)] <- 0
    }

    if(type == "inhomClust"){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)
    }
    # set up the object for bivariate point process
    pp = pp %>%
      mutate(immune = case_when(immune1 == 1 ~ "immune1",
                                immune2 == 1 ~ "immune2",
                                TRUE ~ "background"))




  }else{
    sim_object = CreateSimulationObject(sims = 1, cell_types = 1, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    sim_object = GenerateCellPositivity(sim_object,
                                        k = 25,
                                        sdmin = .7, sdmax = .71,
                                        step_size = 0.1, cores = 1,
                                        probs = c(0.0, 1))

    if(type == "inhomClust"){
      sim_object = GenerateHoles(sim_object, step_size = 0.1, cores = 1)
    }
    # get dataframe with the info you want
    pp = CreateSpatialList(sim_object, single_df = TRUE) %>%
      rename(immune = `Cell 1 Assignment`)

    phat = sum(pp$immune)/nrow(pp)

    if(phat > abundance){
      nhat = nrow(pp)
      nthin = round((phat - abundance) * nhat)

      indices = which(pp$immune == 1)
      pp$immune[sample(indices, nthin)] <- 0
    }

    if(type == "inhomClust"){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)

    }

    pp = pp %>%
      mutate(immune = ifelse(immune == 0, "background", "immune"))

  }# end bivariarte = FALSE

  return(spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune)))

}# end function

