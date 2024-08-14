## function for calculating K in simulated data. Takes in full ppp object
## calculates time for the fperm step.
## can then apply this to a list of ppp objects
get_k = function(ppp_obj, r = c(0, .05, .075,.1, .15, .2)){


  # estimate K using tranlational correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
           r = r,
           correction = c("trans"))

  time_k = toc()

  # calculate Kinhom statistic
  tic()
  kinhom = Kinhom(subset(ppp_obj, marks == "immune"),
                  r = r,
                  correction = c("trans"))
  time_kinhom = toc()


  # calculate fperm statistic
  tic()
  fperm = Kest(ppp_obj,
               r = r,
               correction = c("trans"))
  time_fperm = toc()

  # calculate fperm statistic on 50% thinned data
  tic()
  ppp_obj_thin = rthin(ppp_obj, P = .5)
  fperm_thin = Kest(ppp_obj_thin,
                    r = r,
                    correction = c("trans"))
  time_fperm_thin = toc()



  # calculate perm statistic
  kf = function(obj){
    kdf = Kcross(obj, i = "immune", j = "immune",
           r = r,
           correction = c("trans"))

    as_tibble(kdf) %>% filter(r %in% c(.05,.075, .1)) %>% select(r, trans)
  }


  tic()
  perms = rlabel(ppp_obj, nsim = 10000)
  kperm = map_dfr(perms, kf)
  kperm = kperm %>% group_by(r) %>% summarise(trans = mean(trans)) %>% ungroup()
  time_perm = toc()


  res = tibble(r = c(.05,.075, .1),
               ktheo = filter(as_tibble(k), r %in% c(.05,.075, .1))$theo,
               khat = filter(as_tibble(k),r %in% c(.05,.075, .1))$trans,
               kinhomhat = filter(as_tibble(kinhom), r %in% c(.05,.075, .1))$trans,
               kinhomtheo = filter(as_tibble(kinhom), r %in% c(.05,.075, .1))$theo,
               kfperm = filter(as_tibble(fperm), r %in% c(.05,.075, .1))$trans,
               kperm = kperm$trans,
               kfperm_thin = filter(as_tibble(fperm_thin), r %in% c(.05,.075, .1))$trans,
               time_kinhom = time_kinhom$toc - time_kinhom$tic,
               time_k = time_k$toc - time_k$tic,
               time_fperm = (time_fperm$toc - time_fperm$tic) + time_k,
               time_fperm_thin = (time_fperm_thin$toc - time_fperm_thin$tic) + time_k,
               time_kperm = (time_perm$toc - time_perm$tic) + time_k)


  as.matrix(res)
}

