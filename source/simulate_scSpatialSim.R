

# define function for simulating univariate mIF data. Returns object with and without holes
## lambda_n: intensity for background cells
## lambda_m: intensity for marker positive cells
## holes: should an image be simulated with or without holes
## type: defines the distribution of the point process- homogeneous, inhomogeneous, or clustered
sim_scSpatial <- function(lambda_n,
                          abundance,
                          type = c("homClust", "inhomClust")){


  wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  sim_object = CreateSimulationObject(sims = 1, cell_types = 1, window = wm)
  sim_object = GenerateSpatialPattern(sim_object,
                                      lambda = lambda_n/100)

  if(type == "inhomClust"){
    sim_object = GenerateHoles(sim_object, step_size = 0.1, cores = 1)
  }

  sim_object = GenerateCellPositivity(sim_object,
                                      k = 25,
                                      sdmin = .7, sdmax = .71,
                                      step_size = 0.1, cores = 1,
                                      probs = c(0.0, 1))


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


  return(spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune)))

}

