# ppp_obj: ppp dataset
# rvalue: radius used for Ripley's K
# nperm: number of permutation
get_permutation_distribution = function(ppp_obj, rvalue, bivariate = FALSE) {

  npts = npoints(ppp_obj)
  W = Window(ppp_obj)
  areaW = spatstat.geom:::area(W)

  pp_df = as.data.frame(ppp_obj)
  e = edge.Trans(ppp_obj)

  W = as.matrix(dist(as.matrix(select(pp_df, x, y))))
  W = ifelse(W <= rvalue, 1, 0)
  diag(W) = 0
  Wr = W*e
  R0 = sum(Wr)
  R1 = sum(Wr^2)
  R2 = sum(rowSums(Wr)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  if(bivariate){
    m1 = sum(ppp_obj$marks == "immune1")
    m2 = sum(ppp_obj$marks == "immune2")

    f1 = m1*m2/npts/(npts-1)
    f2 = f1*(m1+m2-2)/(npts-2)
    f3 = f1*(m1-1)*(m2-1)/(npts-2)/(npts-3)

    Kmat = Wr[which(ppp_obj$marks == "immune1"),which(ppp_obj$marks == "immune2")]
    K = areaW*sum(Kmat)/m1/m2 # Ripley's K based on translation correction
    mu_K = areaW*R0/npts/(npts-1) # expectation
    var_K = areaW^2*(R1*f1 + R2*f2 + R3*f3)/m1/m1/m2/m2 - mu_K^2   # variance

  }else{
    m = sum(ppp_obj$marks == "immune")
    npairs =  npts * (npts - 1)
    f1 = m*(m-1)/npts/(npts-1)
    f2 = f1*(m-2)/(npts-2)
    f3 = f2*(m-3)/(npts-3)

    Kmat = Wr[which(ppp_obj$marks == "immune"),which(ppp_obj$marks == "immune")]
    K = areaW*sum(Kmat)/m/(m-1) # Ripley's K based on translation correction
    mu_K = areaW*R0/npts/(npts-1) # expectation
    var_K = areaW^2*(2*R1*f1 + 4*R2*f2 + R3*f3)/m/m/(m-1)/(m-1) - mu_K^2   # variance

  }# end ifelse statement


  Z_k = (K-mu_K) / sqrt(var_K) # Test statistic
  pval_appx = pnorm(-Z_k) # approximated p-value based on normal distribution

  result = tibble(
    r = rvalue,
    khat = K,
    kepd = mu_K, # K expectation under permutation distributions\
    kvpd = var_K,
    Zpd = Z_k, # test statistic
    pvalue = min(1, pval_appx)
  )

  return(result)

}

