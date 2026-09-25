# R port of the local posterior probability method of Pejaver et al. (2022),
# translated from the original MATLAB code:
#   github.com/vpejaver/clingen-svi-comp_calibration @ c56d32d4db45389bbd386c909cd49a625e3e2dc7
#   (directory local_posterior_probability/, MIT license)
# Each function names the .m file it reproduces. The exact semantics are preserved:
#   - window [thr - hw, thr + hw] is inclusive; hw grows by successive addition of 'increment'
#   - edge correction 'c' and the minimum-ClinVar / gnomAD-fraction window conditions
#   - posterior = pos / (pos + w * neg), with w fixed (not recomputed in the bootstrap)
#   - threshold selection and one-sided bound identical to get_all_thresholds.m and
#     get_discounted_thresholds.m
# Implementation difference (not a difference in results): window counts use binary search on
# sorted vectors (findInterval) and are vectorised over all thresholds at once, instead of
# looping over each threshold with find().

# main.m, lines 123-128: posterior probabilities that define each evidence strength.
# Order: 1 = VeryStrong, 2 = Strong, 3 = Moderate, 4 = Supporting.
posterior_thresholds <- function(alpha = 0.0441, c = 1124) {
  j <- 1:4
  lr <- c^(1 / 2^(j - 1))
  list(
    Post_p = lr * alpha / ((lr - 1) * alpha + 1),
    Post_b = lr * (1 - alpha) / ((lr - 1) * (1 - alpha) + 1),
    LR = lr
  )
}

# main.m, line 158: weight of one negative example so that the data match the prior alpha
negative_weight <- function(y, alpha) (1 - alpha) * sum(y == 1) / (sum(y == 0) * alpha)

# main.m, line 161: thresholds = unique values of x and g plus floor(min) and ceil(max), descending
score_grid <- function(x, g) {
  s <- c(x, g)
  sort(unique(c(s, floor(min(s)), ceiling(max(s)))), decreasing = TRUE)
}

# get_both_local_posteriors.m
local_posteriors <- function(x, y, g, thrs, w, minpoints = 100, gft = 0.03, increment = 0.1) {
  xp <- sort(x[y == 1]); xn <- sort(x[y == 0]); gs <- sort(g)
  maxthrs <- max(thrs); minthrs <- min(thrs); lg <- length(g)
  # number of elements of sorted v in [lo, hi]
  cnt <- function(v, lo, hi) findInterval(hi, v) - findInterval(lo, v, left.open = TRUE)

  n <- length(thrs)
  post <- numeric(n); finalwindow <- numeric(n)
  todo <- seq_len(n)
  hw <- 0
  while (length(todo)) {
    t  <- thrs[todo]
    lo <- t - hw; hi <- t + hw
    pos <- cnt(xp, lo, hi); neg <- cnt(xn, lo, hi)
    cc <- rep(1, length(todo))
    up <- hi > maxthrs; dn <- !up & lo < minthrs
    cc[up] <- (maxthrs - lo[up]) / (hi[up] - lo[up])
    cc[dn] <- (hi[dn] - minthrs) / (hi[dn] - lo[dn])
    if (any(cc <= 0)) stop("Problem wih computing c")
    ok <- (pos + neg >= cc * minpoints) & !(cnt(gs, lo, hi) < gft * cc * lg)
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

# get_all_thresholds.m: first index with posterior < Post(j), minus 1.
# 'posteriors' is a matrix (rows = replicates) or a vector.
all_thresholds <- function(posteriors, thrs, Post) {
  if (is.null(dim(posteriors))) posteriors <- matrix(posteriors, nrow = 1)
  out <- matrix(NA_real_, nrow(posteriors), length(Post))
  for (b in seq_len(nrow(posteriors))) {
    for (j in seq_along(Post)) {
      first_below <- which(posteriors[b, ] < Post[j])[1]
      if (!is.na(first_below) && first_below - 1 > 0) out[b, j] <- thrs[first_below - 1]
    }
  }
  out
}

# get_discounted_thresholds.m: one-sided bound (default 5 %) over the bootstrap replicates
discounted_thresholds <- function(thresh, B, discountonesided = 0.05, type = c("pathogenic", "benign")) {
  type <- match.arg(type)
  vapply(seq_len(ncol(thresh)), function(j) {
    invalids <- sum(is.na(thresh[, j]))
    if (invalids > discountonesided * B) return(NA_real_)
    t <- sort(thresh[!is.na(thresh[, j]), j], decreasing = (type == "pathogenic"))
    t[floor(discountonesided * B) - invalids + 1]
  }, numeric(1))
}

# One bootstrap replicate (get_both_bootstrapped_posteriors.m): resample (x, y) jointly, g separately
bootstrap_once <- function(x, y, g, thrs, w, minpoints, gft, increment) {
  qx <- sample.int(length(x), length(x), replace = TRUE)
  qg <- sample.int(length(g), length(g), replace = TRUE)
  local_posteriors(x[qx], y[qx], g[qg], thrs, w, minpoints, gft, increment)$post
}

# Full main.m for one tool. Returns everything needed for tables and figures.
# Reproducibility: the original uses rng('shuffle'); here a fixed L'Ecuyer-CMRG seed makes the
# parallel replicates reproducible.
calibrate <- function(x, y, g, alpha = 0.0441, c = 1124, B = 10000, discountonesided = 0.05,
                      minpoints = 100, gft = 0.03, increment = 0.1,
                      workers = 10, seed = 20260925, pct_band = 5) {
  pt <- posterior_thresholds(alpha, c)
  w <- negative_weight(y, alpha)
  thrs <- score_grid(x, g)

  lp <- local_posteriors(x, y, g, thrs, w, minpoints, gft, increment)
  post_p <- lp$post
  post_b <- 1 - rev(post_p)                       # main.m, line 171

  boot_p <- NULL
  if (B > 0) {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, seed)
    parallel::clusterExport(cl, c("local_posteriors", "bootstrap_once"), envir = environment(calibrate))
    boot_p <- do.call(rbind, parallel::parLapply(cl, seq_len(B), function(i, x, y, g, thrs, w, mp, gf, inc)
      bootstrap_once(x, y, g, thrs, w, mp, gf, inc), x, y, g, thrs, w, minpoints, gft, increment))
  }
  boot_b <- if (!is.null(boot_p)) 1 - boot_p[, ncol(boot_p):1, drop = FALSE] else NULL

  thr_p <- all_thresholds(post_p, thrs, pt$Post_p)[1, ]
  thr_b <- all_thresholds(post_b, rev(thrs), pt$Post_b)[1, ]
  bthr_p <- bthr_b <- disc_p <- disc_b <- NULL
  if (!is.null(boot_p)) {
    bthr_p <- all_thresholds(boot_p, thrs, pt$Post_p)
    bthr_b <- all_thresholds(boot_b, rev(thrs), pt$Post_b)
    disc_p <- discounted_thresholds(bthr_p, B, discountonesided, "pathogenic")
    disc_b <- discounted_thresholds(bthr_b, B, discountonesided, "benign")
  }

  # Curves for the figures (plot_posterior_wrapper.m, shaded = 'CI', onesided = 5):
  # 5th percentile of the bootstrap posteriors; MATLAB prctile = R quantile(type = 5)
  band <- function(m) if (is.null(m)) rep(NA_real_, length(thrs))
                      else apply(m, 2, quantile, probs = pct_band / 100, type = 5, names = FALSE)
  curves <- data.frame(score = thrs, post_p = post_p, post_p_p5 = band(boot_p),
                       post_b = rev(post_b), post_b_p5 = rev(band(boot_b)),
                       halfwindow = lp$finalwindow)

  lev <- c("VeryStrong", "Strong", "Moderate", "Supporting")
  list(
    params = list(alpha = alpha, c = c, B = B, discountonesided = discountonesided,
                  minpoints = minpoints, gft = gft, increment = increment, seed = seed,
                  n_plp = sum(y == 1), n_blb = sum(y == 0), n_gnomad = length(g), w = w),
    levels = data.frame(level = lev, LR = pt$LR, Post_p = pt$Post_p, Post_b = pt$Post_b),
    thresholds = data.frame(level = lev,
                            PP3_point = thr_p, PP3_discounted = if (is.null(disc_p)) NA else disc_p,
                            BP4_point = thr_b, BP4_discounted = if (is.null(disc_b)) NA else disc_b),
    curves = curves,
    boot_thresholds = if (is.null(bthr_p)) NULL else
      setNames(as.data.frame(cbind(bthr_p, bthr_b)), c(paste0("PP3_", lev), paste0("BP4_", lev)))
  )
}
