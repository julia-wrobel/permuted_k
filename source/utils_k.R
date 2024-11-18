## function for calculating K in simulated data. Takes in full ppp object
## calculates time for the fperm step.
## can then apply this to a list of ppp objects
get_k = function(ppp_obj,
                 rvec = c(0, .05, .075,.1, .15, .2),
                 nperm = 1000){


  ################################################################################
  ################################################################################
  # estimate K using tranlational correction
  tic()
  k = Kcross(ppp_obj, i = "immune", j = "immune",
             r = rvec,
             correction = c("trans"))
  time_k = toc()

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
get_k_power = function(ppp_obj, rvec = c(0, 0.25, 0.5, 1)){

  ################################################################################
  ################################################################################
  # https://www.esaim-ps.org/articles/ps/pdf/2013/01/ps120027.pdf
  # paper for asymptotic pvalue for Ripley's K under homogeneity
  # estimate K using translation correction
  tic()
  k = Kest(subset(ppp_obj, marks == "immune"),
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
    mutate(method = "kamp")
  time_kamp = toc()


  ################################################################################
  ################################################################################
  # calculate kamp statistic and variance on 50% litened data
  tic()
  ppp_obj_lite = rthin(ppp_obj, P = .5)
  kamplite = map_dfr(rvec, get_permutation_distribution, ppp_obj = ppp_obj_lite) %>%
    mutate(method = "kamplite")
  time_kamplite = toc()


  times = c((time_kamp$toc - time_kamp$tic),
            (time_kamplite$toc - time_kamplite$tic))


  ################################################################################
  ################################################################################
  # aggregate data
  res = bind_rows(kamp, kamplite) %>%
    mutate(time = rep(times, each = length(rvec)))

  return(res)

}






## this function is for getting the variance and hypothesis test for each statistic
## Not gonna do Kinhom here.
get_k_power_permOnly = function(ppp_obj, rvec = c(0, 0.25, 0.5, 1), nperm = 1000){

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


  res2 = res %>%
    mutate(Z = (khat- expectation) / sqrt(var),
           pvalue = pnorm(-Z),
           method = "kperm approx")

  res = bind_rows(res, res2) %>%
    select(r, khat, expectation, var, Z, pvalue, method, time)


  return(res)

}




get_coxPH = function(num_subj, beta, rhoval, k_df){
  beta_vec <- c(beta, -1, 1)    # Coefficients for covariates
  lambda <- 0.01          # Baseline hazard rate
  censoring_rate <- 0.3   # Proportion of censored observations


  # X1 based on KAMP, so that is the true value
  X1 = k_df$kamp

  cov_matrix = matrix(c(1, rhoval,rhoval, 1), nrow = 2)
  X2 =  rnorm(num_subj)
  X = cbind(X1, X2)
  X2 = (X %*% chol(cov_matrix))[,2]

  cov_matrix = matrix(c(1, 0.2,0.2, 1), nrow = 2)
  X3 =  rnorm(num_subj)
  X = cbind(X1, X3)
  X3 = (X %*% chol(cov_matrix))[,2]


  X <- cbind(X1, X2, X3)
  eta <- X %*% beta_vec

  # Simulate survival times based on the exponential distribution
  # Hazard function: h(t|X) = h0(t) * exp(eta)
  # Survival times are drawn from the survival function
  # S(t|X) = exp(-H0(t) * exp(eta)), where H0(t) is the cumulative baseline hazard
  U <- runif(num_subj)
  T <- -log(U) / (lambda * exp(eta))  # Inverse transform sampling

  # Simulate censoring times
  C <- rexp(num_subj, rate = censoring_rate)

  # Observed times and censoring indicators
  observed_time <- pmin(T, C)
  status <- as.numeric(T <= C)  # 1 if event occurred, 0 if censored

  # Create a data frame
  cox_df <- data.frame(
    time = observed_time,
    status = status,
    kamp = k_df$kamp,
    k = k_df$k,
    X2 = X2
  )

  # run cox models for both k and kamp
  mod_k = coxph(Surv(time, status) ~ k + X2 + X3, data = cox_df)
  mod_kamp = coxph(Surv(time, status) ~ kamp + X2 + X3, data = cox_df)

  # extract beta, se beta, CI, p-value, which model it is from
  bind_rows(tidy(mod_k, conf.int = TRUE) %>% mutate(method = "k"),
            tidy(mod_kamp, conf.int = TRUE) %>% mutate(method = "kamp"))
}




