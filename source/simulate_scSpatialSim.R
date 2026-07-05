

# define function for simulating univariate mIF data. Returns object with and without holes
## lambda_n: intensity for background cells
## lambda_m: intensity for marker positive cells
## holes: should an image be simulated with or without holes
## type: defines the distribution of the point process- homogeneous, inhomogeneous, or clustered
##   "homClust"/"inhomClust": cell positivity is spatially clustered (kernel-based)
##   "inhom": inhomogeneous background (holes) but cell positivity is unclustered (random)
##   "bothClust": background and immune cells are both spatially clustered via
##     independent kernels (same mechanism as immune clustering), no holes
##   "bothClustDep": same as "bothClust" but the immune kernel is dependent on
##     (correlated with) the background kernel instead of independent -
##     clusters overlap/coincide rather than occupying separate locations
sim_scSpatial <- function(lambda_n,
                          abundance,
                          type = c("homClust", "inhomClust", "inhom", "bothClust", "bothClustDep"),
                          bivariate = FALSE){


  wm <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))

  if(bivariate){
    sim_object = CreateSimulationObject(sims = 1, cell_types = 2, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    if(type == "inhom"){
      sim_object = GenerateCellPositivity(sim_object,
                                          sdmin = .7, sdmax = .71,
                                          step_size = 0.1, cores = 1,
                                          probs = c(0.0, 1, 1),
                                          no_kernel = TRUE)
    }else{
      sim_object = GenerateCellPositivity(sim_object,
                                          k = 25,
                                          sdmin = .7, sdmax = .71,
                                          step_size = 0.1, cores = 1,
                                          probs = c(0.0, 1, 1))
    }

    if(type %in% c("inhomClust", "inhom")){
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

    if(type %in% c("inhomClust", "inhom")){
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




  }else if(type %in% c("bothClust", "bothClustDep")){
    # background and immune cells are each driven by their own k=25 kernel.
    # "bothClust": random = TRUE keeps the two kernels independent (separate
    #   cluster locations).
    # "bothClustDep": random = FALSE makes the immune kernel dependent on
    #   (correlated with) the background kernel, so clusters overlap/coincide
    #   instead of occupying separate locations. No holes in either case.
    sim_object = CreateSimulationObject(sims = 1, cell_types = 2, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    sim_object = GenerateCellPositivity(sim_object,
                                        k = 25,
                                        sdmin = .7, sdmax = .71,
                                        step_size = 0.1, cores = 1,
                                        probs = c(0.0, 1),
                                        random = (type == "bothClust"))

    # keep both assignments unresolved (don't randomly break ties) so the
    # background cluster can take priority over the immune cluster below
    pp = CreateSpatialList(sim_object, single_df = TRUE, multihit_action = "keep") %>%
      rename(bg_cluster = `Cell 1 Assignment`,
             immune = `Cell 2 Assignment`)

    # drop points that fall in neither kernel field, so both populations
    # show up as visually distinct clusters rather than one clustered
    # population sitting in a diffuse uniform background
    pp = pp %>%
      filter(bg_cluster == 1 | immune == 1) %>%
      mutate(immune = case_when(bg_cluster == 1 ~ 0,
                                immune == 1 ~ 1,
                                TRUE ~ 0))

    phat = sum(pp$immune)/nrow(pp)

    if(phat > abundance){
      nhat = nrow(pp)
      nthin = round((phat - abundance) * nhat)

      indices = which(pp$immune == 1)
      pp$immune[sample(indices, nthin)] <- 0
    }

    pp = pp %>%
      mutate(immune = ifelse(immune == 0, "background", "immune"))

  }else{
    sim_object = CreateSimulationObject(sims = 1, cell_types = 1, window = wm)
    sim_object = GenerateSpatialPattern(sim_object,
                                        lambda = lambda_n/100)

    if(type == "inhom"){
      sim_object = GenerateCellPositivity(sim_object,
                                          sdmin = .7, sdmax = .71,
                                          step_size = 0.1, cores = 1,
                                          probs = c(0.0, 1),
                                          no_kernel = TRUE)
    }else{
      sim_object = GenerateCellPositivity(sim_object,
                                          k = 25,
                                          sdmin = .7, sdmax = .71,
                                          step_size = 0.1, cores = 1,
                                          probs = c(0.0, 1))
    }

    if(type %in% c("inhomClust", "inhom")){
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

    if(type %in% c("inhomClust", "inhom")){
      pp = pp %>%
        rename(hole = `Hole Assignment`) %>%
        filter(hole == "Keep") %>%
        select(-hole)

    }

    pp = pp %>%
      mutate(immune = ifelse(immune == 0, "background", "immune"))

  }# end bivariate = FALSE

  return(spatstat.geom::ppp(pp$x, pp$y, window = wm,  marks = factor(pp$immune)))

}# end function

