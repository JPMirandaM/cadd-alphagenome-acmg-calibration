# Computationally faster, numerically identical version of the bootstrap in pejaver_port.R.
# The original local_posteriors() sorts the resampled scores in every replicate. Because scores
# take few distinct values, window counts can be obtained from per-value counts (tabulate) and
# cumulative sums over the sorted distinct values v:  #{s <= hi} = cum[findInterval(hi, v)] and
# #{s < lo} = cum[findInterval(lo, v, left.open = TRUE)], which are exactly the counts computed on
# the sorted raw vectors. The resampling itself is unchanged (sample.int on x and g, in the same
# order), so with the same seed every replicate, posterior and threshold is identical to the
# original implementation. Verified in auxiliary/A14 on the Pejaver et al. CADD v1.6 data.
# Requires pejaver_port.R to be sourced first (uses posterior_thresholds, all_thresholds, ...).

# Same algorithm as local_posteriors(), from counts per distinct value
local_posteriors_counts <- function(cp, cn, cg, v, thrs, w, minpoints = 100, gft = 0.03, increment = 0.1, lg) {
  Cp <- c(0, cumsum(cp)); Cn <- c(0, cumsum(cn)); Cg <- c(0, cumsum(cg))
  cnt <- function(C, lo, hi) C[findInterval(hi, v) + 1L] - C[findInterval(lo, v, left.open = TRUE) + 1L]
  maxthrs <- max(thrs); minthrs <- min(thrs)
  n <- length(thrs); post <- numeric(n); finalwindow <- numeric(n)
  todo <- seq_len(n); hw <- 0
  while (length(todo)) {
    t <- thrs[todo]; lo <- t - hw; hi <- t + hw
    pos <- cnt(Cp, lo, hi); neg <- cnt(Cn, lo, hi)
    cc <- rep(1, length(todo))
    up <- hi > maxthrs; dn <- !up & lo < minthrs
    cc[up] <- (maxthrs - lo[up]) / (hi[up] - lo[up])
    cc[dn] <- (hi[dn] - minthrs) / (hi[dn] - lo[dn])
    if (any(cc <= 0)) stop("Problem wih computing c")
    ok <- (pos + neg >= cc * minpoints) & !(cnt(Cg, lo, hi) < gft * cc * lg)
    if (any(ok)) {
      idx <- todo[ok]
      post[idx] <- pos[ok] / (pos[ok] + w * neg[ok])
      finalwindow[idx] <- hw
      todo <- todo[!ok]
    }
    hw <- hw + increment
  }
  list(post = post, finalwindow = finalwindow)
}

# Integer codes of x and g on the common sorted distinct values
prepare_codes <- function(x, g) {
  v <- sort(unique(c(x, g)))
  list(v = v, xc = match(x, v), gc = match(g, v))
}

# One variant-level bootstrap replicate: identical draws to bootstrap_once()
bootstrap_once_fast <- function(xc, y, gc, v, thrs, w, minpoints, gft, increment) {
  nx <- length(xc); ng <- length(gc); k <- length(v)
  qx <- sample.int(nx, nx, replace = TRUE)
  qg <- sample.int(ng, ng, replace = TRUE)
  xq <- xc[qx]; yq <- y[qx]
  local_posteriors_counts(tabulate(xq[yq == 1], k), tabulate(xq[yq == 0], k), tabulate(gc[qg], k),
                          v, thrs, w, minpoints, gft, increment, lg = ng)$post
}

# calibrate() of pejaver_port.R using the fast replicate functions. Same parallel scheme
# (parLapply over 1..B with clusterSetRNGStream(seed)), so with the same seed and number of workers
# the replicates are identical to calibrate(). Optional gene-level bootstrap (gx, gg = gene labels).
calibrate_fast <- function(x, y, g, alpha = 0.0441, c = 1124, B = 10000, discountonesided = 0.05,
                           minpoints = 100, gft = 0.03, increment = 0.1, workers = 10, seed = 20260925,
                           pct_band = 5, cluster = FALSE, gx = NULL, gg = NULL) {
  pt <- posterior_thresholds(alpha, c); w <- negative_weight(y, alpha); thrs <- score_grid(x, g)
  lp <- local_posteriors(x, y, g, thrs, w, minpoints, gft, increment)          # point estimate: original code
  post_p <- lp$post; post_b <- 1 - rev(post_p)
  pc <- prepare_codes(x, g)
  gx_idx <- if (cluster) unname(split(seq_along(gx), gx)) else NULL
  gg_idx <- if (cluster) unname(split(seq_along(gg), gg)) else NULL
  cl <- parallel::makeCluster(workers); on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, seed)
  parallel::clusterExport(cl, c("local_posteriors_counts", "bootstrap_once_fast", "bootstrap_once_cluster_fast"),
                          envir = globalenv())
  boot_p <- do.call(rbind, parallel::parLapply(cl, seq_len(B), function(i, xc, y, gc, v, thrs, w, mp, gf, inc, cluster, gxi, ggi)
    if (cluster) bootstrap_once_cluster_fast(xc, y, gc, v, thrs, w, mp, gf, inc, gxi, ggi)
    else bootstrap_once_fast(xc, y, gc, v, thrs, w, mp, gf, inc),
    pc$xc, y, pc$gc, pc$v, thrs, w, minpoints, gft, increment, cluster, gx_idx, gg_idx))
  boot_b <- 1 - boot_p[, ncol(boot_p):1, drop = FALSE]
  bthr_p <- all_thresholds(boot_p, thrs, pt$Post_p); bthr_b <- all_thresholds(boot_b, rev(thrs), pt$Post_b)
  band <- function(m) apply(m, 2, quantile, probs = pct_band / 100, type = 5, names = FALSE)
  lev <- c("VeryStrong", "Strong", "Moderate", "Supporting")
  list(
    levels = data.frame(level = lev, LR = pt$LR, Post_p = pt$Post_p, Post_b = pt$Post_b),
    thresholds = data.frame(level = lev,
      PP3_point = all_thresholds(post_p, thrs, pt$Post_p)[1, ],
      PP3_discounted = discounted_thresholds(bthr_p, B, discountonesided, "pathogenic"),
      BP4_point = all_thresholds(post_b, rev(thrs), pt$Post_b)[1, ],
      BP4_discounted = discounted_thresholds(bthr_b, B, discountonesided, "benign"),
      PP3_boot_NA_frac = colMeans(is.na(bthr_p)), BP4_boot_NA_frac = colMeans(is.na(bthr_b))),
    curves = data.frame(score = thrs, post_p = post_p, post_p_p5 = band(boot_p), post_b = rev(post_b),
                        post_b_p5 = rev(band(boot_b)), halfwindow = lp$finalwindow),
    boot_thresholds = setNames(as.data.frame(cbind(bthr_p, bthr_b)), c(paste0("PP3_", lev), paste0("BP4_", lev))),
    params = list(alpha = alpha, c = c, B = B, seed = seed, workers = workers, cluster_bootstrap = cluster,
                  minpoints = minpoints, gft = gft, increment = increment, n_plp = sum(y == 1),
                  n_blb = sum(y == 0), n_gnomad = length(g), w = w, implementation = "pejaver_port_fast.R"))
}

# One gene-level (cluster) bootstrap replicate (sensitivity analysis S3)
bootstrap_once_cluster_fast <- function(xc, y, gc, v, thrs, w, minpoints, gft, increment, gx_idx, gg_idx) {
  k <- length(v)
  sx <- sample.int(length(gx_idx), length(gx_idx), replace = TRUE)
  sg <- sample.int(length(gg_idx), length(gg_idx), replace = TRUE)
  ix <- unlist(gx_idx[sx], use.names = FALSE); ig <- unlist(gg_idx[sg], use.names = FALSE)
  xq <- xc[ix]; yq <- y[ix]
  local_posteriors_counts(tabulate(xq[yq == 1], k), tabulate(xq[yq == 0], k), tabulate(gc[ig], k),
                          v, thrs, w, minpoints, gft, increment, lg = length(ig))$post
}
