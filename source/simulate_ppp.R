# define function for simulating univariate mIF data. Returns object with and without holes
## lambda_n: intensity for background cells
## lambda_nm: intensity for marker positive cells
## holes: should an image be simulated with or without holes
## type: defines the distribution of the point process- homogeneous, inhomogeneous, or clustered
mxsim_univariate <- function(lambda_n,
                             lambda_nm, # needs to be divisible by 5
                             type = c("hom", "inhom", "homClust", "inhomClust")){



  wm <- circle_window()

  if(type %in% c("inhom", "inhomClust")){
    lams <- list(function(x,y){
      #lambda_nm*exp(-.6*x^2 + 0.5*y)
      lambda_nm*x^2
    }
    # log quadratic trend
    ,
    function(x,y){
      #lambda_n*exp(-0.6*x^2+0.5*y)
      lambda_n*x^2
    }
    # log linear trend
    )
  }


  if(type == "hom"){
    # homogeneous background and immune
    pp_obj = rmpoispp(c(lambda_nm, lambda_n), types = c("immune", "tumor"),
                      win = wm)
  }else if(type == "inhom"){
    # inhomogeneous background and inhomogeneous immune
    pp_obj = rmpoispp(lams, types = c("immune", "tumor"),
                      win = wm)

  }else if(type == "homClust"){
    # homogeneous background, clustered immune
    pp_obj_tumor = rpoispp(lambda_n, win = wm)
    marks(pp_obj_tumor) = "tumor"
    pp_obj_clust = rMatClust(5, 0.05, lambda_nm / 5, win = wm)
    marks(pp_obj_clust) = "immune"

    pp_df = bind_rows(as_tibble(pp_obj_tumor), as_tibble(pp_obj_clust)) %>%
      mutate(marks = factor(marks))

    pp_obj = ppp(pp_df$x, pp_df$y,wm, marks = pp_df$marks)


  }else if(type == "inhomClust"){
    # inhomogeneous background, clustered immune
    pp_obj_tumor = rpoispp(lams[[2]], win = wm)
    marks(pp_obj_tumor) = "tumor"
    pp_obj_clust = rMatClust(5, 0.05, lambda_nm / 5, win = wm)
    marks(pp_obj_clust) = "immune"
    pp_df = bind_rows(as_tibble(pp_obj_tumor), as_tibble(pp_obj_clust)) %>%
      mutate(marks = factor(marks))

    pp_obj = ppp(pp_df$x, pp_df$y,wm, marks = pp_df$marks)
  }


  # define holes here
  pp_obj_holes = simulate_holes(pp_obj)
  list(full = pp_obj, holes = pp_obj_holes)
}



### simulate holes
simulate_holes <- function(ppp_obj){
  #create random number of points
  verts = round(rnorm(1, mean = 25, sd = 8))
  #create verts number of scaling factors
  scalex = abs(rnorm(verts, 5, 2))
  scaley = abs(rnorm(verts, 5, 2))
  #create verts number of translational factors
  transx = runif(1, -1, 1)
  transy = runif(1, -1, 1)
  #finally, generate vert number of x and y and perform scaling
  x = runif(verts, -1, 1)/scalex + transx
  y = runif(verts, -1, 1)/scaley + transy
  #if x and y are outside our window, limit them to window's edge
  x[x > 1] = 1; x[x < -1] = -1
  y[y > 1] = 1; y[y < -1] = -1
  #create a convex hull using our new x and y points
  tmp_w = convexhull.xy(x=x, y=y)
  #create a new window with our hole in it, and convert to a mask of TRUE and FALSE
  w2 = owin(c(-1,1), c(-1,1), poly = list(x=tmp_w$bdry[[1]]$x, y=tmp_w$bdry[[1]]$y))
  tmp = as.mask(w2, dimyx = c(100, 100))
  # 100 x 100 image, all TRUE
  #subtract our original TMA core window and make sure it contains logical values for a mask
  w <- owin(c(-1,1), c(-1,1), mask=matrix(TRUE, 100,100))
  X <- raster.x(w)
  Y <- raster.y(w)
  wm <- owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))

  blh = tmp$m - t(as.matrix(wm))
  blh = sapply(as.data.frame(blh), as.logical)
  #create new TMA core window for this new core and assign the subtracted window to it
  wmp = owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))
  wmp$m = blh
  #assign new window, then reapply the original window to remove cells that are in the hole
  Window(ppp_obj) <- wmp
  Window(ppp_obj) <- wm
  return(ppp_obj)
}

circle_window <- function(){
  w <- owin(c(-1,1), c(-1,1), mask=matrix(TRUE, 100,100))
  X <- raster.x(w)
  Y <- raster.y(w)
  wm <- owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))
  wm
}
