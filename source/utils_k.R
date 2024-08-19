## function for calculating K in simulated data. Takes in full ppp object
## calculates time for the fperm step.
## can then apply this to a list of ppp objects
get_k = function(ppp_obj, rvec = c(0, .05, .075,.1, .15, .2), nperm = 10000){


  ################################################################################
  ################################################################################
  # estimate K using tranlational correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
           r = rvec,
           correction = c("trans"))

  time_k = toc()

  ################################################################################
  ################################################################################
  # calculate Kinhom statistic
  tic()
  kinhom = Kinhom(subset(ppp_obj, marks == "immune"),
                  r = rvec,
                  correction = c("trans"))
  time_kinhom = toc()


  ################################################################################
  ################################################################################
  # calculate kepd statistic
  tic()
  kepd = Kest(ppp_obj,
               r = rvec,
               correction = c("trans"))
  time_kepd = toc()

  ################################################################################
  ################################################################################
  # calculate kepd statistic on 50% thinned data
  tic()
  ppp_obj_thin = rthin(ppp_obj, P = .5)
  kepd_thin = Kest(ppp_obj_thin,
                    r = rvec,
                    correction = c("trans"))
  time_kepdThin = toc()

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
  kperm = kperm %>% group_by(r) %>% summarise(trans = mean(trans)) %>% ungroup()
  time_perm = toc()


  ################################################################################
  ################################################################################
  res = tibble(r = rvec,
               ktheo = filter(as_tibble(k), r %in% rvec)$theo,
               khat = filter(as_tibble(k),r %in% rvec)$trans,
               kinhomhat = filter(as_tibble(kinhom), r %in% rvec)$trans,
               kinhomtheo = filter(as_tibble(kinhom), r %in% rvec)$theo,
               kepd = filter(as_tibble(kepd), r %in% rvec)$trans,
               kperm = kperm$trans,
               kepdThin = filter(as_tibble(kepd_thin), r %in% rvec)$trans,
               time_kinhom = time_kinhom$toc - time_kinhom$tic,
               time_k = time_k$toc - time_k$tic,
               time_kepd = (time_kepd$toc - time_kepd$tic) + time_k,
               time_kepdThin = (time_kepdThin$toc - time_kepdThin$tic) + time_k,
               time_kperm = (time_perm$toc - time_perm$tic) + time_k)


  as.matrix(res)
}



## this function is for getting the variance and hypothesis test for each statistic
## Not gonna do Kinhom here.
get_k_power = function(ppp_obj, rvec = c(0, .15), nperm = 10000){

  ################################################################################
  ################################################################################
  # estimate K using translation correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
             r = rvec,
             correction = c("trans"),
             var.approx = TRUE)

  # using envelope approach to get pvalue
  simulated_k <- envelope(subset(ppp_obj, marks == "immune"), Kest, nsim = nperm,
                          correction = c("trans"))

  as_tibble(simulated_k) %>%
    filter(r %in% rvec)
    group_by(r) %>%
    mutate()
  pvalue <- sum(k$trans > simulated_k$obs) / nperm
  variance_CSR <- k$var
  time_k = toc()


  ################################################################################
  ################################################################################
  # calculate kepd statistic and variance
  tic()
  kepd = map_dfr(rvec, get_permutation_distribution, ppp_obj = ppp_obj, variance = TRUE) %>%
    select(r, var = kvpd, pvalue, expectation = kepd)
  time_kepd = toc()


  ################################################################################
  ################################################################################
  # calculate kepd statistic and variance on 50% thinned data
  tic()
  ppp_obj_thin = rthin(ppp_obj, P = .5)
  kepdThin = map_dfr(rvec, get_permutation_distribution, ppp_obj = ppp_obj_thin, variance = TRUE) %>%
    select(r, var = kvpd, pvalue, expectation = kepd)
  time_kepdThin = toc()


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
                                              pvalue = sum((khat-trans) >= 0)/nperm,
                                              expectation = mean(trans)) %>% ungroup()
  time_perm = toc()

  ################################################################################
  ################################################################################
  # aggregate data
  res = tibble(r = rvec,
               ktheo = filter(as_tibble(k), r %in% rvec)$theo,
               khat = filter(as_tibble(k),r %in% rvec)$trans,
               kinhomhat = filter(as_tibble(kinhom), r %in% rvec)$trans,
               kinhomtheo = filter(as_tibble(kinhom), r %in% rvec)$theo,
               kepd = filter(as_tibble(kepd), r %in% rvec)$trans,
               kperm = kperm$trans,
               kepdThin = filter(as_tibble(kepd_thin), r %in% rvec)$trans,
               time_kinhom = time_kinhom$toc - time_kinhom$tic,
               time_k = time_k$toc - time_k$tic,
               time_kepd = (time_kepd$toc - time_kepd$tic) + time_k,
               time_kepdThin = (time_kepdThin$toc - time_kepdThin$tic) + time_k,
               time_kperm = (time_perm$toc - time_perm$tic) + time_k)


  as.matrix(res)
}






