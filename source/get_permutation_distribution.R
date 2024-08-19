# ppp_obj: ppp dataset
# rvalue: radius used for Ripley's K
# nperm: number of permutation
get_permutation_distribution = function(ppp_obj, rvalue, variance = TRUE) {
  if(!variance){
    kepd = Kest(ppp_obj,
                r = rvalue,
                correction = c("trans"))
    return(kepd)
  }else{
    m = sum(ppp_obj$marks == "immune")
    npts = npoints(ppp_obj)
    npairs =  npts * (npts - 1)
    W = Window(ppp_obj)
    areaW = area(W)

    rmax = diameter(W) + 1
    close = closepairs(ppp_obj, rmax, what="all")
    DIJ = close$d

    DIJ.mat = matrix(0, npts, npts)
    DIJ.mat[lower.tri(DIJ.mat)] = DIJ[1:(npairs/2)]
    DIJ.mat[upper.tri(DIJ.mat)] = t(DIJ.mat)[upper.tri(DIJ.mat)]

    gW = setcov(W) # what is this
    edgewt = edge.Trans(dx=close$dx, dy=close$dy, W=W, paired=TRUE, gW=gW, give.rmax=TRUE)

    trans.mat = matrix(0, npts, npts)
    trans.mat[lower.tri(trans.mat)] = edgewt[1:(npairs/2)]
    trans.mat[upper.tri(trans.mat)] = t(trans.mat)[upper.tri(trans.mat)]

    D = matrix(as.numeric(DIJ.mat<=rvalue), npts, npts)
    Wr = D*trans.mat

    R0 = sum(Wr)
    R1 = sum(Wr^2)
    R2 = sum(rowSums(Wr)^2) - R1
    R3 = R0^2 - 2*R1 - 4*R2

    f1 = m*(m-1)/npts/(npts-1)
    f2 = f1*(m-2)/(npts-2)
    f3 = f2*(m-3)/(npts-3)

    indices = which(ppp_obj$marks=="immune")
    Kmat = Wr[indices,indices]

    khat = areaW*sum(Kmat)/m/(m-1) # Ripley's K based on translation correction

    mu_k = areaW*R0/npts/(npts-1) # expectation
    var_k = areaW^2*(2*R1*f1 + 4*R2*f2 + R3*f3)/m/m/(m-1)/(m-1) - mu_k^2   # variance

    Z_k = (khat-mu_k) / sqrt(var_k) # Test statistic
    pval_appx = pnorm(-Z_k) # approximated p-value based on normal distribution


    result = tibble(
      r = rvalue,
      khat = khat,
      kepd = mu_k, # K expectation under permutation distributions\
      kvpd = var_k,
      Zpd = Z_k, # test statistic
      pvalue = min(1, pval_appx)
    )

    return(result)
  }

}

