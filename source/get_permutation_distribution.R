# ppp_obj: ppp dataset
# rvalue: radius used for Ripley's K
# bivariate: logical, whether to use bivariate K
#
# Uses closepairs() (KD-tree) to find only the O(n*k) pairs within radius r,
# avoiding the full O(n^2) distance matrix. Edge corrections are computed
# without ever building an n x n matrix:
#   - rectangular windows: vectorized formula per pair
#   - polygonal/mask windows: single 2D FFT autocorrelation of the mask (same
#     pixel approximation spatstat uses internally in edge.Trans), then sparse
#     lookup by pair displacement
# Produces the same E(K) and Var(K) as the original dense implementation.
get_permutation_distribution = function(ppp_obj, rvalue, bivariate = FALSE) {

  npts  = npoints(ppp_obj)
  W_win = Window(ppp_obj)
  areaW = spatstat.geom:::area(W_win)

  # KD-tree range query: all ordered pairs (i,j), i!=j, d(i,j) <= rvalue
  cp    = closepairs(ppp_obj, rmax = rvalue, what = "all", twice = TRUE)
  i_idx = cp$i
  j_idx = cp$j

  # --- edge corrections for pairs within radius r only ---
  if (W_win$type == "rectangle") {
    W_width  = diff(W_win$xrange)
    W_height = diff(W_win$yrange)
    area_int = pmax(0, W_width - abs(cp$dx)) * pmax(0, W_height - abs(cp$dy))
    e_vals   = ifelse(area_int > 0, areaW / area_int, 0)

  } else {
    # Convert polygon to mask (matches what edge.Trans does internally).
    # Compute 2D autocorrelation via separable 1D FFTs once, then look up
    # each pair's displacement — never materialise an n x n matrix.
    W_mask = as.mask(W_win)
    m_num  = as.numeric(W_mask$m)          # 0/1 vector (row-major after matrix())
    ny     = nrow(W_mask$m); nx = ncol(W_mask$m)
    xstep  = W_mask$xstep;  ystep = W_mask$ystep

    # 2D DFT: FFT each column, then each row
    F_m = apply(matrix(m_num, ny, nx), 2, fft)
    F_m = t(apply(F_m, 1, fft))

    # 2D circular autocorrelation scaled to spatial area units:
    # ac[m1+1, m2+1] = area(W ∩ shift(W, -m2*xstep, -m1*ystep))
    ac_step = apply(Mod(F_m)^2, 2, fft, inverse = TRUE) / ny
    ac      = Re(t(apply(ac_step, 1, fft, inverse = TRUE))) / nx * xstep * ystep

    # For pair (i,j) with displacement (dx, dy) we want
    # area(W ∩ shift(W, +dx, +dy)) = ac[((-di) %% ny)+1, ((-dj) %% nx)+1]
    di   = round(cp$dy / ystep)
    dj   = round(cp$dx / xstep)
    ridx = ((-di) %% ny) + 1L
    cidx = ((-dj) %% nx) + 1L

    overlap = ac[cbind(ridx, cidx)]
    e_vals  = ifelse(overlap > 0, areaW / overlap, 0)
  }

  # --- R statistics via sparse row sums, no n x n matrix built ---
  R0 = sum(e_vals)
  R1 = sum(e_vals^2)

  row_sums = numeric(npts)
  if (length(i_idx) > 0) {
    agg = tapply(e_vals, i_idx, sum)
    row_sums[as.integer(names(agg))] = agg
  }
  R2 = sum(row_sums^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  marks_vec = as.character(ppp_obj$marks)

  if (bivariate) {
    m1    = sum(marks_vec == "immune1")
    m2    = sum(marks_vec == "immune2")
    f1    = m1*m2/npts/(npts-1)
    f2    = f1*(m1+m2-2)/(npts-2)
    f3    = f1*(m1-1)*(m2-1)/(npts-2)/(npts-3)
    keep  = marks_vec[i_idx] == "immune1" & marks_vec[j_idx] == "immune2"
    K     = areaW * sum(e_vals[keep]) / m1 / m2
    mu_K  = areaW * R0 / npts / (npts-1)
    var_K = areaW^2 * (R1*f1 + R2*f2 + R3*f3) / m1/m1/m2/m2 - mu_K^2
  } else {
    m     = sum(marks_vec == "immune")
    f1    = m*(m-1)/npts/(npts-1)
    f2    = f1*(m-2)/(npts-2)
    f3    = f2*(m-3)/(npts-3)
    keep  = marks_vec[i_idx] == "immune" & marks_vec[j_idx] == "immune"
    K     = areaW * sum(e_vals[keep]) / m / (m-1)
    mu_K  = areaW * R0 / npts / (npts-1)
    var_K = areaW^2 * (2*R1*f1 + 4*R2*f2 + R3*f3) / m/m/(m-1)/(m-1) - mu_K^2
  }

  Z_k       = (K - mu_K) / sqrt(var_K)
  pval_appx = pnorm(-Z_k)

  tibble(
    r           = rvalue,
    khat        = K,
    expectation = mu_K,
    var         = var_K,
    Z           = Z_k,
    pvalue      = min(1, pval_appx)
  )
}
