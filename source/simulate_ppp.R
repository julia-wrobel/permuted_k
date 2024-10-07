# define function for simulating univariate mIF data. Returns object with and without holes
## lambda_n: intensity for background cells
## lambda_m: intensity for marker positive cells
## holes: should an image be simulated with or without holes
## type: defines the distribution of the point process- homogeneous, inhomogeneous, or clustered
mxsim_univariate <- function(lambda_n,
                             abundance, # needs to be divisible by 5
                             type = c("hom", "inhom")){



  wm <- spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1))

  lambda_immune = round((lambda_n * abundance)/(1-abundance))
  lambda_background = lambda_n - lambda_immune

  if(type %in% c("inhom")){
    lams <- list(function(x,y){
      lambda_immune*5*x^2
    }
    ,
    function(x,y){
      lambda_background*5*x^2
    }
    # log linear trend
    )
  }


  if(type == "hom"){
    # homogeneous background and immune
    pp_obj = rmpoispp(c(lambda_immune, lambda_background), types = c("immune", "background"),
                      win = wm)
  }else if(type == "inhom"){
    # inhomogeneous background and inhomogeneous immune
    pp_obj = rmpoispp(lams, types = c("immune", "background"),
                      win = wm)

  }


  return(pp_obj)
}


