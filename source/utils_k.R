## function for calculating K in simulated data. Takes in full ppp object
## calculates time for the fperm step.
## can then apply this to a list of ppp objects
get_k = function(ppp_obj,
                 #rvec = c(0, .05, .075,.1, .15, .2),
                 nperm = 1000){


  ################################################################################
  ################################################################################
  # estimate K using tranlational correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
             #r = rvec,
             correction = c("trans"))
  time_k = toc()

  # pull vector of r values
  rvec = k$r

  k = k %>%
    as_tibble() %>%
    mutate(method = "k") %>%
    select(r, csr = theo, trans, method)



  ################################################################################
  ################################################################################
  # calculate Kinhom statistic
  tic()
  kinhom = Kinhom(subset(ppp_obj, marks == "immune"),
                  r = rvec,
                  correction = c("trans")) %>%
    as_tibble() %>% mutate(method = "kinhom") %>%
    select(r, csr = theo, trans, method)
  time_kinhom = toc()


  ################################################################################
  ################################################################################
  # calculate KAMP statistic
  tic()
  kamp = Kest(ppp_obj,
              r = rvec,
              correction = c("trans")) %>%
    as_tibble() %>%
    mutate(method = "kamp",
           csr = trans,
           trans = k$trans) %>%
    select(r, csr, trans, method)
  time_kamp = toc()


  ################################################################################
  ################################################################################
  # calculate kamp statistic on 50% thinned data
  tic()
  ppp_obj_lite = rthin(ppp_obj, P = .5)
  kamp_lite = Kest(ppp_obj_lite,
                   r = rvec,
                   correction = c("trans")) %>%
    as_tibble() %>%
    mutate(method = "kamp_lite",
           csr = trans,
           trans = k$trans) %>%
    select(r, csr, trans, method)
  time_kamplite = toc()

  ################################################################################
  ################################################################################
  # calculate perm statistic
  kf = function(obj){
    kdf = Kcross(obj, i = "immune", j = "immune",
                 r = rvec,
                 correction = c("trans"))

    as_tibble(kdf) %>% filter(r %in% rvec) %>% select(r, trans)
  }


  tic()
  perms = rlabel(ppp_obj, nsim = nperm)
  kperm = map_dfr(perms, kf)
  kperm = kperm %>% group_by(r) %>% summarise(csr = mean(trans)) %>% ungroup() %>%
    mutate(trans = k$trans,
           method = "perm") %>%
    select(r, csr, trans, method)
  time_perm = toc()


  ################################################################################
  ################################################################################

  times = c(time_k$toc - time_k$tic,
             time_kinhom$toc - time_kinhom$tic,
             time_kamp$toc - time_kamp$tic,
             time_kamplite$toc - time_kamplite$tic,
             time_perm$toc - time_perm$tic)


  ################################################################################
  ################################################################################
  # aggregate data
  res = bind_rows(k, kinhom, kamp, kamp_lite, kperm) %>%
    mutate(time = rep(times, each = length(rvec)))

  return(res)



}


## this function is for getting the variance and hypothesis test for each statistic
## Not gonna do Kinhom here.
get_k_power = function(ppp_obj, rvec = c(0, .15)){

  ################################################################################
  ################################################################################
  # https://www.esaim-ps.org/articles/ps/pdf/2013/01/ps120027.pdf
  # paper for asymptotic pvalue for Ripley's K under homogeneity
  # estimate K using translation correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
             r = rvec,
             correction = c("trans"))
  time_k = toc()
  # Not using anyliteg for hypothesis testing for K under theoretical CSR. Need to build this in later.
  # get variance based on block bootstrap for use in confidence intervals
  # not the same as the permutation variance, let' s
  # I don't litek this is right
  #var_k = envelope(ppp_obj, fun = Kcross, i = "immune",
   #        j = "immune", r = rvec, correction = c("trans"), global = FALSE,
    #       nsim = 999, alternative = "greater")

  ################################################################################
  ################################################################################
  # calculate kamp statistic and variance
  tic()
  kamp = map_dfr(rvec, get_permutation_distribution, ppp_obj = ppp_obj) %>%
    select(r, var = kvpd, pvalue, expectation = kamp) %>%
    mutate(method = "kamp")
  time_kamp = toc()


  ################################################################################
  ################################################################################
  # calculate kamp statistic and variance on 50% litened data
  tic()
  ppp_obj_lite = rthin(ppp_obj, P = .5)
  kamplite = map_dfr(rvec, get_permutation_distribution, ppp_obj = ppp_obj_lite) %>%
    select(r, var = kvpd, pvalue, expectation = kamp) %>%
    mutate(method = "kamplite")
  time_kamplite = toc()


  times = c((time_kamp$toc - time_kamp$tic),
            (time_kamplite$toc - time_kamplite$tic))


  ################################################################################
  ################################################################################
  # aggregate data
  res = bind_rows(kamp, kamplite) %>%
    mutate(khat = rep(k$trans, times = 2),
           time = rep(times, each = length(rvec)))

  return(res)

}






## this function is for getting the variance and hypothesis test for each statistic
## Not gonna do Kinhom here.
get_k_power_permOnly = function(ppp_obj, rvec = c(0, .15), nperm = 10000){

  ################################################################################
  ################################################################################
  # https://www.esaim-ps.org/articles/ps/pdf/2013/01/ps120027.pdf
  # paper for asymptotic pvalue for Ripley's K under homogeneity
  # estimate K using translation correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
             r = rvec,
             correction = c("trans"))
  time_k = toc()
  # Not using anyliteg for hypothesis testing for K under theoretical CSR. Need to build this in later.
  # get variance based on block bootstrap for use in confidence intervals
  # not the same as the permutation variance, let' s
  # I don't litek this is right
  #var_k = envelope(ppp_obj, fun = Kcross, i = "immune",
  #        j = "immune", r = rvec, correction = c("trans"), global = FALSE,
  #       nsim = 999, alternative = "greater")



  ################################################################################
  ################################################################################
  # calculate perm statistic
  kf = function(obj){
    kdf = Kcross(obj, i = "immune", j = "immune",
                 r = rvec,
                 correction = c("trans"))

    as_tibble(kdf) %>% filter(r %in% rvec) %>% select(r, trans) %>%
      mutate(khat = k$trans)
  }


  tic()
  perms = rlabel(ppp_obj, nsim = nperm)
  kperm = map_dfr(perms, kf)
  kperm = kperm %>% group_by(r) %>% summarise(var = var(trans),
                                              pvalue = sum(trans >= khat)/nperm,
                                              expectation = mean(trans)) %>% ungroup() %>%
    mutate(method = "kperm")
  time_perm = toc()


  times = c((time_perm$toc - time_perm$tic))


  ################################################################################
  ################################################################################
  # aggregate data
  res = kperm %>%
    mutate(khat = rep(k$trans, times = 1),
           time = rep(times, each = length(rvec)))

  return(res)

}






