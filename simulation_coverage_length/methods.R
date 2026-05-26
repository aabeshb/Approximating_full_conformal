
# ============================================================
# - beta_true = sqrt(10) * u, where u is uniform on unit sphere in R^p
# - ols_beta uses Moore–Penrose inverse when system is undetermined / rank-deficient
# ============================================================

# ---------- Helper: union length from intervals ----------
interval_union_length <- function(intervals) {
  if (length(intervals) == 0) return(0)
  M <- do.call(rbind, intervals)
  M <- M[order(M[, 1], M[, 2]), , drop = FALSE]
  total <- 0
  curL <- M[1, 1]; curR <- M[1, 2]
  for (i in seq_len(nrow(M))[-1]) {
    L <- M[i, 1]; R <- M[i, 2]
    if (L <= curR) {
      curR <- max(curR, R)
    } else {
      total <- total + (curR - curL)
      curL <- L; curR <- R
    }
  }
  total + (curR - curL)
}

# ---------- Random beta: sqrt(10) * (uniform unit vector) ----------
# Uniform on the unit sphere via z ~ N(0, I_p), u = z / ||z||
rand_beta_unit_scaled <- function(p, scale = sqrt(10)) {
  z <- rnorm(p)
  nz <- sqrt(sum(z^2))
  scale * (z / nz)
}

# ---------- OLS helpers (NO intercept), MP inverse fallback ----------
ols_beta <- function(X, y) {
  # If system is underdetermined or rank-deficient, use Moore–Penrose inverse
  qrx <- qr(X)
  if (nrow(X) < ncol(X) || qrx$rank < ncol(X)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Please install MASS for ginv(): install.packages('MASS')")
    }
    return(drop(MASS::ginv(X) %*% y))
  }
  # Otherwise, use fast QR-based least squares
  drop(qr.coef(qrx, y))
}

ols_predict <- function(X, beta) as.numeric(X %*% beta)

# ---------- Grid construction ----------
make_grid <- function(y, M, pad = 0.02) {
  a <- min(y)
  b <- max(y)
  if (!is.finite(a) || !is.finite(b) || a == b) {
    a <- min(y); b <- max(y)
    if (a == b) { a <- a - 1; b <- b + 1 }
  }
  rng <- b - a
  a <- a - pad * rng
  b <- b + pad * rng
  seq(a, b, length.out = M)
}

######################### Approximate methods ##############################

# --------------------------- Rounding: approximate  ---------------------------
gridcp_length <- function(X, y, x_new, grid, alpha) {
  n <- nrow(X); M <- length(grid)
  
  # Voronoi cells
  mids <- (grid[-M] + grid[-1]) / 2
  lefts  <- c(-Inf, mids)
  rights <- c(mids,  Inf)
  
  intervals <- list()
  
  # robust quantile probability
  prob <- min(1, (1 - alpha) * (1 + 1/n))
  
  for (m in seq_len(M)) {
    X_aug <- rbind(X, matrix(x_new, nrow = 1))
    y_aug <- c(y, grid[m])
    beta  <- ols_beta(X_aug, y_aug)
    
    if (anyNA(beta)) next
    
    mu_m  <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta))
    resid <- abs(y - ols_predict(X, beta))
    if (anyNA(resid)) next
    
    q_m <- as.numeric(quantile(resid, probs = prob, type = 1, names = FALSE))
    
    L <- max(mu_m - q_m, lefts[m])
    R <- min(mu_m + q_m, rights[m])
    if (L < R) intervals[[length(intervals) + 1L]] <- c(L, R)
  }
  
  if (length(intervals) == 0) return(0)
  interval_union_length(intervals)
}

# ---------------------------- MEMBERSHIP -----------------------------------
gridcp_contains <- function(X, y, x_new, y_true, grid, alpha) {
  n <- nrow(X)
  prob <- min(1, (1 - alpha) * (1 + 1/n))
  
  m <- which.min(abs(grid - y_true))
  
  X_aug <- rbind(X, matrix(x_new, nrow = 1))
  y_aug <- c(y, grid[m])
  beta  <- ols_beta(X_aug, y_aug)
  
  if (anyNA(beta)) return(NA)
  
  mu_m  <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta))
  resid <- abs(y - ols_predict(X, beta))
  if (anyNA(resid)) return(NA)
  
  q_m <- as.numeric(quantile(resid, probs = prob, type = 1, names = FALSE))
  abs(y_true - mu_m) <= q_m
}

# --------------------- Deletion: approximate  -----------------------------------
naive_cp_interval_linear_abs <- function(X, y, x_new, alpha) {
  n <- nrow(X)
  
  beta_full <- ols_beta(X, y)
  if (anyNA(beta_full)) return(c(NA_real_, NA_real_))
  
  mu_new <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta_full))
  q_hat = sort(abs(y - ols_predict(X,beta_full)))[ceiling((1-alpha)*(n+1))]
  
  L <- mu_new - q_hat
  U <- mu_new + q_hat
  
  return(c(L, U))
}

naive_cp_length_linear_abs <- function(X, y, x_new, alpha) {
  int <- naive_cp_interval_linear_abs(X, y, x_new, alpha)
  if (anyNA(int)) return(NA_real_)
  max(0, int[2] - int[1])
}

naive_cp_contains_linear_abs <- function(X, y, x_new, y_true, alpha) {
  int <- naive_cp_interval_linear_abs(X, y, x_new, alpha)
  if (anyNA(int)) return(NA)
  (y_true >= int[1]) && (y_true <= int[2])
}


# --------------------- One-step update: approximate ---------------------------
onegrad_set_linear_abs <- function(X, y, x_new, alpha, eta = 0.1, eps = 1e-12) {
  n <- nrow(X)
  
  beta_hat <- ols_beta(X, y)
  if (anyNA(beta_hat)) {
    return(list(intervals = list(), beta_hat = NA))
  }
  
  d    <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta_hat))
  yhat <- ols_predict(X, beta_hat)
  e    <- y - yhat
  
  kappa     <- as.numeric((eta / (n + 1)) * (X %*% x_new))
  kappa_new <- as.numeric((eta / (n + 1)) * sum(x_new^2))
  ccoef     <- abs(1 - kappa_new)  # slope of test score S_new(y)=ccoef*|y-d|
  
  # If c ~ 0, then S_new(y) ~ 0 for all y, so the set is essentially all R
  if (!is.finite(ccoef) || ccoef < eps) {
    return(list(intervals = list(c(-Inf, Inf)),
                beta_hat = beta_hat, d = d, c = ccoef,
                kappa = kappa, kappa_new = kappa_new, e = e))
  }
  
  # quantile index matching quantile(..., type=1)
  prob <- min(1, (1 - alpha) * (1 + 1/n))
  k <- ceiling(prob * n)
  k <- max(1L, min(n, k))
  
  # 2n breakpoints from the two linear equations:
  #   e_i - kappa_i t = + c t  => t = e_i/(kappa_i + c)
  #   e_i - kappa_i t = - c t  => t = e_i/(kappa_i - c)
  den_plus  <- kappa + ccoef
  den_minus <- kappa - ccoef
  
  t_plus  <- ifelse(abs(den_plus)  > eps, e / den_plus,  NA_real_)
  t_minus <- ifelse(abs(den_minus) > eps, e / den_minus, NA_real_)
  
  B <- c(d, d + t_plus, d + t_minus)
  B <- B[is.finite(B)]
  B <- sort(unique(B))
  
  edges <- c(-Inf, B, Inf)
  
  intervals <- list()
  for (j in seq_len(length(edges) - 1L)) {
    L <- edges[j]
    R <- edges[j + 1L]
    
    # pick a representative point strictly inside (L, R)
    if (is.infinite(L) && is.finite(R)) {
      mid <- R - 1
    } else if (is.finite(L) && is.infinite(R)) {
      mid <- L + 1
    } else if (is.finite(L) && is.finite(R)) {
      mid <- L + (R - L) / 2
    } else {
      mid <- 0
    }
    
    t <- mid - d
    s_new <- ccoef * abs(t)
    s_i   <- abs(e - kappa * t)
    
    count_lt <- sum(s_i < s_new)
    if (count_lt <= (k - 1L)) {
      intervals[[length(intervals) + 1L]] <- c(L, R)
    }
  }
  
  list(intervals = intervals,
       beta_hat = beta_hat, d = d, c = ccoef,
       kappa = kappa, kappa_new = kappa_new, e = e,
       k = k, prob = prob)
}

onegrad_length_linear_abs <- function(X, y, x_new, alpha, eta = 0.1) {
  obj <- onegrad_set_linear_abs(X, y, x_new, alpha, eta = eta)
  if (length(obj$intervals) == 0) return(0)
  interval_union_length(obj$intervals)
}

onegrad_contains_linear_abs <- function(X, y, x_new, y_true, alpha, eta = 0.1) {
  obj <- onegrad_set_linear_abs(X, y, x_new, alpha, eta = eta)
  ints <- obj$intervals
  if (length(ints) == 0) return(FALSE)
  for (int in ints) {
    if (y_true >= int[1] && y_true <= int[2]) return(TRUE)
  }
  FALSE
}



################# Tournament corrected methods ############################
round_one_length <- function(X, y, x_new, grid, alpha) {
  n <- nrow(X)
  idx_nearest <- function(v) which.min(abs(grid - v))
  mids <- (grid[-length(grid)] + grid[-1]) / 2
  left_edges <- c(-Inf, mids)
  right_edges <- c(mids, Inf)
  
  tau <- (1 - alpha) * (n + 1)
  K <- n - floor(tau)
  
  total_len <- 0
  
  for (m in seq_along(grid)) {
    y_m <- grid[m]
    evx <- numeric(0); evd <- integer(0)
    
    for (i in 1:n) {
      y_mod <- y
      y_mod[i] <- grid[idx_nearest(y[i])]
      
      X_aug <- rbind(X, matrix(x_new, nrow = 1))
      y_aug <- c(y_mod, y_m)
      beta <- ols_beta(X_aug, y_aug)
      
      mu_i <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta))
      r_i  <- abs(y[i] - ols_predict(matrix(X[i,], nrow = 1), beta))
      L <- mu_i - r_i; R <- mu_i + r_i
      if (L < R) { evx <- c(evx, L, R); evd <- c(evd, +1L, -1L) }
    }
    
    if (length(evx) == 0) next
    events_df <- aggregate(delta ~ x, data = data.frame(x = evx, delta = evd), sum)
    events_df <- events_df[order(events_df$x), , drop = FALSE]
    
    coverage <- 0L
    prev_x <- NULL
    segs <- list()
    
    for (k in seq_len(nrow(events_df))) {
      xk <- events_df$x[k]
      if (!is.null(prev_x) && coverage >= K && xk > prev_x) {
        segs[[length(segs) + 1]] <- c(prev_x, xk)
      }
      coverage <- coverage + events_df$delta[k]
      prev_x <- xk
    }
    
    if (length(segs) > 0) {
      for (s in segs) {
        Lc <- max(s[1], left_edges[m]); Rc <- min(s[2], right_edges[m])
        if (Lc < Rc) total_len <- total_len + (Rc - Lc)
      }
    }
  }
  total_len
}

# ---------- "Rounding: corrected membership ----------
round_one_contains <- function(X, y, x_new, y_true, grid, alpha) {
  n <- nrow(X)
  idx_nearest <- function(v) which.min(abs(grid - v))
  m <- idx_nearest(y_true)
  
  tau <- (1 - alpha) * (n + 1)
  out_count <- 0L
  
  for (i in 1:n) {
    y_mod <- y
    y_mod[i] <- grid[idx_nearest(y[i])]
    X_aug <- rbind(X, matrix(x_new, nrow = 1))
    y_aug <- c(y_mod, grid[m])
    beta <- ols_beta(X_aug, y_aug)
    
    mu_i <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta))
    r_i  <- abs(y[i] - ols_predict(matrix(X[i,], nrow = 1), beta))
    out_count <- out_count + as.integer(abs(y_true - mu_i) > r_i)
    if (out_count >= tau) return(FALSE)
  }
  TRUE
}

######################### Leave-one-out-cross-conformal ###################

# ---------- helper: merge closed intervals ----------
merge_closed_intervals <- function(intervals) {
  if (length(intervals) == 0) return(list())
  
  M <- do.call(rbind, intervals)
  M <- M[order(M[, 1], M[, 2]), , drop = FALSE]
  
  out <- list()
  curL <- M[1, 1]
  curR <- M[1, 2]
  
  if (nrow(M) >= 2) {
    for (i in 2:nrow(M)) {
      L <- M[i, 1]
      R <- M[i, 2]
      
      # closed intervals: merge if they overlap or touch
      if (L <= curR) {
        curR <- max(curR, R)
      } else {
        out[[length(out) + 1L]] <- c(curL, curR)
        curL <- L
        curR <- R
      }
    }
  }
  
  out[[length(out) + 1L]] <- c(curL, curR)
  out
}


# ---------- precompute the n leave-one-out intervals ----------
loo_cross_ols_precompute <- function(X, y, x_new) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  x_new <- as.numeric(x_new)
  
  n <- nrow(X)
  p <- ncol(X)
  
  stopifnot(length(y) == n, length(x_new) == p)
  
  mu_new_loo <- numeric(n)
  R <- numeric(n)
  lo <- numeric(n)
  hi <- numeric(n)
  
  for (i in seq_len(n)) {
    beta_minus <- ols_beta(X[-i, , drop = FALSE], y[-i])
    
    mu_new <- ols_predict(matrix(x_new, nrow = 1), beta_minus)
    mu_i   <- ols_predict(matrix(X[i, ], nrow = 1), beta_minus)
    
    R_i <- abs(y[i] - mu_i)
    
    mu_new_loo[i] <- mu_new
    R[i] <- R_i
    lo[i] <- mu_new - R_i
    hi[i] <- mu_new + R_i
  }
  
  list(
    lo = lo,
    hi = hi,
    mu_new_loo = mu_new_loo,
    R = R
  )
}


# ---------- exact LOO cross-conformal set ----------
# returns a list of closed intervals
loo_cross_conformal_set_linear_abs <- function(X, y, x_new, alpha) {
  stopifnot(0 <= alpha, alpha <= 1)
  
  tmp <- loo_cross_ols_precompute(X, y, x_new)
  lo <- tmp$lo
  hi <- tmp$hi
  n  <- length(lo)
  
  # Need count >= tau, where (1 + count)/(n + 1) > alpha
  tau <- floor(alpha * (n + 1))
  
  if (tau <= 0) return(list(c(-Inf, Inf)))
  if (tau > n)  return(list())
  
  u <- sort(unique(c(lo, hi)))
  
  starts <- tabulate(match(lo, u), nbins = length(u))
  ends   <- tabulate(match(hi, u), nbins = length(u))
  
  # number of covering intervals on (u_j, u_{j+1})
  count_after_u <- cumsum(starts - ends)
  
  pieces <- list()
  
  # coverage exactly at endpoint u_j
  count_before_u <- c(0L, head(count_after_u, -1L))
  depth_at_u <- count_before_u + starts
  
  idx_pt <- which(depth_at_u >= tau)
  if (length(idx_pt) > 0L) {
    pieces <- c(pieces, lapply(idx_pt, function(j) c(u[j], u[j])))
  }
  
  # coverage on open intervals (u_j, u_{j+1})
  if (length(u) >= 2L) {
    idx_seg <- which(head(count_after_u, -1L) >= tau)
    if (length(idx_seg) > 0L) {
      pieces <- c(pieces, lapply(idx_seg, function(j) c(u[j], u[j + 1L])))
    }
  }
  
  merge_closed_intervals(pieces)
}


# ---------- containment ----------
loo_cross_conformal_contains_linear_abs <- function(X, y, x_new, y_true, alpha) {
  stopifnot(0 <= alpha, alpha <= 1)
  
  tmp <- loo_cross_ols_precompute(X, y, x_new)
  lo <- tmp$lo
  hi <- tmp$hi
  n  <- length(lo)
  
  tau <- floor(alpha * (n + 1))
  if (tau <= 0) return(TRUE)
  if (tau > n)  return(FALSE)
  
  sum(lo <= y_true & y_true <= hi) >= tau
}


# ---------- total length ----------
loo_cross_conformal_length_linear_abs <- function(X, y, x_new, alpha) {
  stopifnot(0 <= alpha, alpha <= 1)
  
  tmp <- loo_cross_ols_precompute(X, y, x_new)
  lo <- tmp$lo
  hi <- tmp$hi
  n  <- length(lo)
  
  tau <- floor(alpha * (n + 1))
  if (tau <= 0) return(Inf)
  if (tau > n)  return(0)
  
  u <- sort(unique(c(lo, hi)))
  if (length(u) <= 1L) return(0)
  
  starts <- tabulate(match(lo, u), nbins = length(u))
  ends   <- tabulate(match(hi, u), nbins = length(u))
  
  # depth on each open interval (u_j, u_{j+1})
  depth_open <- head(cumsum(starts - ends), -1L)
  
  sum((u[-1L] - u[-length(u)]) * (depth_open >= tau))
}
######################### Corrected onegrad ################################

corrected_onegrad_object_linear_abs <- function(X, y, x_new, alpha, eta = 0.1, eps = 1e-12) {
  n <- nrow(X)
  p <- ncol(X)
  
  a <- numeric(n)   # a_i = eta/(n+1) * <X_i, x_new>
  g <- numeric(n)   # g_i = eta/(n+1) * ||X_i||^2
  r <- numeric(n)   # r_i^- = Y_i - X_i^T beta_{-i}
  d <- numeric(n)   # d_i   = x_new^T beta_{-i}
  
  b <- as.numeric((eta / (n + 1)) * sum(x_new^2))   # common across i
  breaks <- numeric(0)
  
  for (i in seq_len(n)) {
    beta_minus <- ols_beta(X[-i, , drop = FALSE], y[-i])
    if (anyNA(beta_minus)) {
      return(list(ok = FALSE, intervals = list()))
    }
    
    m_i <- as.numeric(ols_predict(matrix(X[i, ], nrow = 1), beta_minus))
    r[i] <- y[i] - m_i
    d[i] <- as.numeric(ols_predict(matrix(x_new, nrow = 1), beta_minus))
    
    a[i] <- as.numeric((eta / (n + 1)) * sum(X[i, ] * x_new))
    g[i] <- as.numeric((eta / (n + 1)) * sum(X[i, ]^2))
    
    # Breakpoints from
    # |(1-b)(y-d_i) - a_i r_i| = |(1-g_i) r_i - a_i (y-d_i)|
      # Let t = y - d_i. Then:
      # |(1-b)t - a_i r_i| = |(1-g_i)r_i - a_i t|
          #
          # Case 1: (1-b)t - a_i r_i =  (1-g_i)r_i - a_i t
          #         => (1-b+a_i)t = (1-g_i+a_i)r_i
          #
          # Case 2: (1-b)t - a_i r_i = -(1-g_i)r_i + a_i t
          #         => (1-b-a_i)t = (a_i + g_i - 1)r_i
          
          den1 <- 1 - b + a[i]
          den2 <- 1 - b - a[i]
          
          if (abs(den1) > eps) {
            breaks <- c(breaks, d[i] + ((1 - g[i]) + a[i]) * r[i] / den1)
          }
          if (abs(den2) > eps) {
            breaks <- c(breaks, d[i] + (a[i] + g[i] - 1) * r[i] / den2)
          }
  }
  
  B <- sort(unique(breaks[is.finite(breaks)]))
  edges <- c(-Inf, B, Inf)
  
  thresh <- (1 - alpha) * (n + 1)
  intervals <- list()
  
  for (j in seq_len(length(edges) - 1L)) {
    L <- edges[j]
    R <- edges[j + 1L]
    
    # representative point inside (L, R)
    if (is.infinite(L) && is.finite(R)) {
      y_mid <- R - 1
    } else if (is.finite(L) && is.infinite(R)) {
      y_mid <- L + 1
    } else if (is.finite(L) && is.finite(R)) {
      y_mid <- (L + R) / 2
    } else {
      y_mid <- 0
    }
    
    t_vec <- y_mid - d
    
    test_scores  <- abs((1 - b) * t_vec - a * r)
    calib_scores <- abs((1 - g) * r - a * t_vec)
    
    count_bad <- sum(test_scores > calib_scores)
    
    if (count_bad < thresh) {
      intervals[[length(intervals) + 1L]] <- c(L, R)
    }
  }
  
  list(
    ok = TRUE,
    intervals = intervals,
    a = a,
    b = b,
    g = g,
    r = r,
    d = d,
    thresh = thresh
  )
}

corrected_onegrad_length_from_object <- function(obj) {
  if (is.null(obj$ok) || !isTRUE(obj$ok)) return(NA_real_)
  if (length(obj$intervals) == 0) return(0)
  interval_union_length(obj$intervals)
}

corrected_onegrad_contains_from_object <- function(obj, y_true) {
  if (is.null(obj$ok) || !isTRUE(obj$ok)) return(NA)
  
  t_vec <- y_true - obj$d
  
  test_scores  <- abs((1 - obj$b) * t_vec - obj$a * obj$r)
  calib_scores <- abs((1 - obj$g) * obj$r - obj$a * t_vec)
  
  sum(test_scores > calib_scores) < obj$thresh
}

corrected_onegrad_length_linear_abs <- function(X, y, x_new, alpha, eta = 0.1, eps = 1e-12) {
  obj <- corrected_onegrad_object_linear_abs(X, y, x_new, alpha, eta = eta, eps = eps)
  corrected_onegrad_length_from_object(obj)
}

corrected_onegrad_contains_linear_abs <- function(X, y, x_new, y_true, alpha, eta = 0.1, eps = 1e-12) {
  obj <- corrected_onegrad_object_linear_abs(X, y, x_new, alpha, eta = eta, eps = eps)
  corrected_onegrad_contains_from_object(obj, y_true)
}

# ==============================================================================
# Bayesian AOI approximation and tournament version
# ==============================================================================

logsumexp <- function(v) {
  v <- as.numeric(v)
  vmax <- max(v)
  if (!is.finite(vmax)) return(vmax)
  vmax + log(sum(exp(v - vmax)))
}

bayesian_lr <- function(mu_0, Sigma_0, X, y, T = 100, sig_lik_model = 1.0) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install MASS: install.packages('MASS')")
  }

  sigma2 <- sig_lik_model^2
  Sigma0_inv <- solve(Sigma_0)

  Sigma_n <- solve(crossprod(X) / sigma2 + Sigma0_inv)
  mu_n <- drop(Sigma_n %*% (crossprod(X, y) / sigma2 + Sigma0_inv %*% mu_0))

  theta_T <- MASS::mvrnorm(n = T, mu = mu_n, Sigma = Sigma_n)
  if (is.null(dim(theta_T))) {
    theta_T <- matrix(theta_T, ncol = length(mu_0))
  }
  theta_T
}

theta_loo <- function(mu_0, Sigma_0, X, y, T = 100, sig_lik_model = 1.0) {
  n <- length(y)
  d <- length(mu_0)
  theta_nT <- array(0, dim = c(n, T, d))

  for (j in seq_len(n)) {
    keep <- seq_len(n) != j
    theta_nT[j, , ] <- bayesian_lr(
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      X = X[keep, , drop = FALSE],
      y = y[keep],
      T = T,
      sig_lik_model = sig_lik_model
    )
  }

  theta_nT
}

loglik_bayes <- function(y, x, theta, sig_lik_model = 1.0) {
  theta <- as.matrix(theta)
  x <- as.numeric(x)

  mu <- drop(theta %*% x)
  sigma2 <- sig_lik_model^2
  -0.5 * log(2.0 * pi * sigma2) - (y - mu)^2 / (2.0 * sigma2)
}

# Approximate Bayesian AOI bad-event count using samples from pi(theta | D_n).
bayes_aoi_bad_count <- function(X, x_new, y, y_cand, theta_T, sig_lik_model = 1.0) {
  n <- length(y)
  bad_count <- 0L

  logw <- loglik_bayes(y_cand, x_new, theta_T, sig_lik_model = sig_lik_model)
  log_ppd_test_num <- logsumexp(2.0 * logw)

  for (i in seq_len(n)) {
    log_ppd_train_num <- logsumexp(
      logw + loglik_bayes(y[i], X[i, ], theta_T, sig_lik_model = sig_lik_model)
    )
    bad_count <- bad_count + as.integer(log_ppd_test_num < log_ppd_train_num)
  }

  bad_count
}

# Corrected Bayesian tournament bad-event count using samples from pi(theta | D_n \ {Z_j}).
bayes_corrected_bad_count <- function(X, x_new, y, y_cand, theta_nT, sig_lik_model = 1.0) {
  n <- length(y)
  bad_count <- 0L

  for (j in seq_len(n)) {
    theta_j <- theta_nT[j, , ]
    if (is.null(dim(theta_j))) {
      theta_j <- matrix(theta_j, ncol = ncol(X))
    }

    loglik_j <- loglik_bayes(y[j], X[j, ], theta_j, sig_lik_model = sig_lik_model)
    loglik_test <- loglik_bayes(y_cand, x_new, theta_j, sig_lik_model = sig_lik_model)

    log_ppd_test_num <- logsumexp(loglik_j + 2.0 * loglik_test)
    log_ppd_train_num <- logsumexp(2.0 * loglik_j + loglik_test)

    bad_count <- bad_count + as.integer(log_ppd_test_num < log_ppd_train_num)
  }

  bad_count
}

satisfies_quantile_condition <- function(bad_count, n, alpha) {
  bad_count < (1.0 - alpha) * (n + 1)
}

# Grid-based interval approximation for Bayesian set length.
# This is only used for estimating average length in the simulation.
bayes_CI_shrink <- function(X, x_new, y, theta_obj, bad_count_fn,
                            alpha = 0.10,
                            min_y = NULL,
                            max_y = NULL,
                            size = 0.1,
                            shrink = 10,
                            sig_lik_model = 1.0) {
  if (is.null(min_y) || is.null(max_y)) {
    y_min <- min(y)
    y_max <- max(y)
    y_range <- y_max - y_min
    min_y <- y_min - 0.02 * y_range
    max_y <- y_max + 0.02 * y_range
  }

  delta <- (max_y - min_y) / shrink
  left <- min_y
  right <- max_y

  n <- length(y)
  q <- (1.0 - alpha) * (n + 1)
  check <- 0L

  eval_bad <- function(y_val) {
    bad_count_fn(
      X,
      x_new,
      y,
      y_val,
      theta_obj,
      sig_lik_model = sig_lik_model
    )
  }

  if (delta > size) {
    grid <- seq(min_y, max_y + delta, by = delta)

    for (i in seq_along(grid)) {
      y_val <- grid[i]
      bad_count <- eval_bad(y_val)

      if ((bad_count < q) && (check == 0L)) {
        left <- grid[max(i - 1L, 1L)]
        check <- 1L
      }

      if ((y_val > left) && (bad_count >= q) && (check == 1L)) {
        right <- y_val
        break
      }
    }

    delta <- delta / shrink

    while (delta > size) {
      if (check == 1L) {
        grid <- seq(0.0, (shrink + 1) * delta, by = delta)
        left_check <- 0L
        right_check <- 0L

        for (i in seq_along(grid)) {
          y_val <- grid[i]

          if (left_check == 0L) {
            bad_left <- eval_bad(left + y_val)
            if (bad_left < q) {
              left <- left + grid[max(i - 1L, 1L)]
              left_check <- 1L
            }
          }

          if (right_check == 0L) {
            bad_right <- eval_bad(right - y_val)
            if (bad_right < q) {
              right <- right - grid[max(i - 1L, 1L)]
              right_check <- 1L
            }
          }

          if ((left_check == 1L) && (right_check == 1L)) break
        }
      }

      if (check == 0L) {
        grid <- seq(min_y, max_y + delta, by = delta)
        for (i in seq_along(grid)) {
          y_val <- grid[i]
          bad_count <- eval_bad(y_val)

          if ((bad_count < q) && (check == 0L)) {
            left <- grid[max(i - 1L, 1L)]
            check <- 1L
          }

          if ((y_val > left) && (bad_count >= q) && (check == 1L)) {
            right <- y_val
            break
          }
        }
      }

      delta <- delta / shrink
    }

    if (check == 1L) {
      grid <- seq(0.0, (shrink + 1) * size, by = size)
      left_check <- 0L
      right_check <- 0L

      for (i in seq_along(grid)) {
        y_val <- grid[i]

        if (left_check == 0L) {
          bad_left <- eval_bad(left + y_val)
          if (bad_left < q) {
            left <- left + grid[max(i - 1L, 1L)]
            left_check <- 1L
          }
        }

        if (right_check == 0L) {
          bad_right <- eval_bad(right - y_val)
          if (bad_right < q) {
            right <- right - grid[max(i - 1L, 1L)]
            right_check <- 1L
          }
        }

        if ((left_check == 1L) && (right_check == 1L)) break
      }
    }

    if (check == 0L) {
      grid <- seq(min_y, max_y + size, by = size)
      for (i in seq_along(grid)) {
        y_val <- grid[i]
        bad_count <- eval_bad(y_val)

        if ((bad_count < q) && (check == 0L)) {
          left <- grid[max(i - 1L, 1L)]
          check <- 1L
        }

        if ((y_val > left) && (bad_count >= q) && (check == 1L)) {
          right <- y_val
          break
        }
      }
    }
  } else {
    grid <- seq(min_y, max_y + size, by = size)
    check <- 0L

    for (i in seq_along(grid)) {
      y_val <- grid[i]
      bad_count <- eval_bad(y_val)

      if ((bad_count < q) && (check == 0L)) {
        left <- grid[max(i - 1L, 1L)]
        check <- 1L
      }

      if ((y_val > left) && (bad_count >= q) && (check == 1L)) {
        right <- y_val
        break
      }
    }
  }

  if (check == 0L) {
    left <- 0.0
    right <- 0.0
  }

  c(left, right)
}

bayesian_is_results_for_dataset <- function(X, y, x_new, y_true,
                                            alpha = 0.10,
                                            T = 100,
                                            mu_prior = 10.0,
                                            sig_prior = 1.0,
                                            sig_lik_model = 1.0,
                                            min_y_aoi = NULL,
                                            max_y_aoi = NULL,
                                            min_y_corrected = NULL,
                                            max_y_corrected = NULL,
                                            length_grid_size = 0.1,
                                            shrink = 10) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  x_new <- as.numeric(x_new)

  y_min <- min(y)
  y_max <- max(y)
  y_range <- y_max - y_min
  min_y_aoi <- y_min - 0.02 * y_range
  max_y_aoi <- y_max + 0.02 * y_range
  min_y_corrected <- min_y_aoi
  max_y_corrected <- max_y_aoi

  d <- ncol(X)
  n <- nrow(X)

  mu_0 <- rep(mu_prior, d)
  Sigma_0 <- diag(sig_prior^2, d)

  theta_T <- bayesian_lr(
    mu_0 = mu_0,
    Sigma_0 = Sigma_0,
    X = X,
    y = y,
    T = T,
    sig_lik_model = sig_lik_model
  )

  theta_nT <- theta_loo(
    mu_0 = mu_0,
    Sigma_0 = Sigma_0,
    X = X,
    y = y,
    T = T,
    sig_lik_model = sig_lik_model
  )

  aoi_bad <- bayes_aoi_bad_count(
    X = X,
    x_new = x_new,
    y = y,
    y_cand = y_true,
    theta_T = theta_T,
    sig_lik_model = sig_lik_model
  )

  corrected_bad <- bayes_corrected_bad_count(
    X = X,
    x_new = x_new,
    y = y,
    y_cand = y_true,
    theta_nT = theta_nT,
    sig_lik_model = sig_lik_model
  )

  CI_aoi <- bayes_CI_shrink(
    X = X,
    x_new = x_new,
    y = y,
    theta_obj = theta_T,
    bad_count_fn = bayes_aoi_bad_count,
    alpha = alpha,
    min_y = min_y_aoi,
    max_y = max_y_aoi,
    size = length_grid_size,
    shrink = shrink,
    sig_lik_model = sig_lik_model
  )

  CI_corrected <- bayes_CI_shrink(
    X = X,
    x_new = x_new,
    y = y,
    theta_obj = theta_nT,
    bad_count_fn = bayes_corrected_bad_count,
    alpha = alpha,
    min_y = min_y_corrected,
    max_y = max_y_corrected,
    size = length_grid_size,
    shrink = shrink,
    sig_lik_model = sig_lik_model
  )

  data.frame(
    method = c("Bayesian:approx", "Bayesian:corrected"),
    covered = c(
      satisfies_quantile_condition(aoi_bad, n, alpha),
      satisfies_quantile_condition(corrected_bad, n, alpha)
    ),
    length = c(
      CI_aoi[2] - CI_aoi[1],
      CI_corrected[2] - CI_corrected[1]
    )
  )
}
