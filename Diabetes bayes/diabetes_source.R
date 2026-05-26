# ============================================================
# Fong-style Bayesian-lasso conformal prediction in R
# Source file: functions only; no experiment is run when sourced.
#
# Implements only the Bayesian methods:
#   Bayesian:approx     = add-one-in importance-sampling conformal Bayes
#   Bayesian:corrected  = tournament correction
#
# Model:
#   Y_i | X_i, beta, intercept, sigma ~ N(X_i^T beta + intercept, sigma^2)
#   beta_j | b                         ~ Laplace(0, b)
#   b                                  ~ Gamma(1, 1)
#   intercept                          ~ weak N(0,10) prior; data are standardized
#   sigma                              ~ HalfNormal(sigma_prior_scale)
#
# ============================================================

# ------------------------- setup -----------------------------

check_rstan <- function() {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Please install rstan first. See https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started")
  }
}

logsumexp <- function(v) {
  v <- as.numeric(v)
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

col_logsumexp <- function(M) {
  M <- as.matrix(M)
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colLogSumExps(M))
  }
  m <- apply(M, 2, max)
  ans <- m + log(colSums(exp(sweep(M, 2, m, "-"))))
  ans[!is.finite(m)] <- m[!is.finite(m)]
  ans
}


read_rds_if_valid <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  tryCatch(
    readRDS(file),
    error = function(e) {
      warning("Could not read cached file ", file, "; refitting. Error: ", conditionMessage(e))
      NULL
    }
  )
}

save_rds_atomic <- function(object, file) {
  if (is.null(file)) return(invisible(FALSE))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(file, ".tmp.", Sys.getpid())
  saveRDS(object, tmp)
  ok <- file.rename(tmp, file)
  if (!ok) {
    # Fallback for filesystems where file.rename across temp paths can fail.
    saveRDS(object, file)
    if (file.exists(tmp)) unlink(tmp)
  }
  invisible(TRUE)
}

standardize_split <- function(X_train, y_train, X_test, y_test = NULL) {
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)
  if (!is.null(y_test)) y_test <- as.numeric(y_test)

  x_mean <- colMeans(X_train)
  x_sd <- apply(X_train, 2, stats::sd)
  x_sd[!is.finite(x_sd) | x_sd == 0] <- 1

  y_mean <- mean(y_train)
  y_sd <- stats::sd(y_train)
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1

  list(
    X_train = sweep(sweep(X_train, 2, x_mean, "-"), 2, x_sd, "/"),
    y_train = (y_train - y_mean) / y_sd,
    X_test = sweep(sweep(X_test, 2, x_mean, "-"), 2, x_sd, "/"),
    y_test = if (is.null(y_test)) NULL else (y_test - y_mean) / y_sd,
    x_mean = x_mean,
    x_sd = x_sd,
    y_mean = y_mean,
    y_sd = y_sd
  )
}

# ------------------------- Stan model ------------------------

fong_blasso_stan_code <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n, p] X;
  vector[n] y;
  real<lower=0> sigma_prior_scale;
}
parameters {
  vector[p] beta;
  real intercept;
  real<lower=1e-6> b;
  real<lower=1e-6> sigma;
}
model {
  b ~ gamma(1, 1);
  beta ~ double_exponential(0, b);
  intercept ~ normal(0, 10);
  sigma ~ normal(0, sigma_prior_scale);
  y ~ normal(intercept + X * beta, sigma);
}
"

compile_fong_blasso_model <- function() {
  check_rstan()
  rstan::rstan_options(auto_write = TRUE)
  rstan::stan_model(model_code = fong_blasso_stan_code)
}

fit_fong_blasso_mcmc <- function(
    X,
    y,
    stan_model = NULL,
    iter = 3000,
    warmup = 1000,
    chains = 4,
    seed = 1,
    sigma_prior_scale = 1,
    adapt_delta = 0.8,
    max_treedepth = 10,
    refresh = 0,
    keep_fit = FALSE
) {
  check_rstan()
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(y) == n)

  if (is.null(stan_model)) stan_model <- compile_fong_blasso_model()

  init_fun <- function(chain_id = 1) {
    list(
      beta = rep(0, p),
      intercept = mean(y),
      b = 1,
      sigma = max(stats::sd(y), 0.5)
    )
  }

  fit <- tryCatch(
    rstan::sampling(
      object = stan_model,
      data = list(n = n, p = p, X = X, y = y, sigma_prior_scale = sigma_prior_scale),
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed,
      init = init_fun,
      refresh = refresh,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    ),
    error = function(e) {
      stop(
        "Stan sampling failed before returning draws. Re-run with chains=1, refresh=1, ",
        "and check the Stan messages. Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  post <- tryCatch(
    rstan::extract(fit, pars = c("beta", "intercept", "sigma", "b"), permuted = TRUE),
    error = function(e) {
      stop(
        "Stan returned a fit object with no usable posterior samples. ",
        "This usually means all chains failed. Re-run with chains=1 and refresh=1. ",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (is.null(post$beta) || length(post$beta) == 0) {
    stop("No posterior beta draws were returned by Stan.", call. = FALSE)
  }
  beta_draws <- post$beta
  if (is.null(dim(beta_draws))) beta_draws <- matrix(beta_draws, ncol = p)

  out <- list(
    beta = beta_draws,
    beta0 = as.numeric(post$intercept),
    sigma = as.numeric(post$sigma),
    b = as.numeric(post$b),
    iter = iter,
    warmup = warmup,
    chains = chains
  )
  if (keep_fit) out$fit <- fit
  out
}

fit_loo_fong_blasso_mcmc <- function(
    X,
    y,
    stan_model = NULL,
    iter = 3000,
    warmup = 1000,
    chains = 4,
    seed = 1,
    sigma_prior_scale = 1,
    adapt_delta = 0.8,
    max_treedepth = 10,
    refresh = 0,
    cache_dir = NULL,
    verbose = FALSE
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(stan_model)) stan_model <- compile_fong_blasso_model()
  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  out <- vector("list", n)
  for (j in seq_len(n)) {
    cache_file <- if (is.null(cache_dir)) NULL else file.path(cache_dir, sprintf("loo_%03d.rds", j))
    cached <- read_rds_if_valid(cache_file)
    if (!is.null(cached)) {
      out[[j]] <- cached
      if (verbose && (j == 1L || j == n || j %% 25L == 0L)) {
        cat(sprintf("Loaded cached LOO posterior %d/%d\n", j, n))
      }
      next
    }

    if (verbose && (j == 1L || j == n || j %% 25L == 0L)) {
      cat(sprintf("Fitting LOO posterior %d/%d\n", j, n))
    }
    keep <- seq_len(n) != j
    out[[j]] <- fit_fong_blasso_mcmc(
      X = X[keep, , drop = FALSE],
      y = y[keep],
      stan_model = stan_model,
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed + j,
      sigma_prior_scale = sigma_prior_scale,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh,
      keep_fit = FALSE
    )
    save_rds_atomic(out[[j]], cache_file)
  }
  out
}

# ---------------- likelihood evaluations ---------------------

loglik_blasso <- function(y_val, x, samples) {
  x <- as.numeric(x)
  beta <- as.matrix(samples$beta)
  mu <- samples$beta0 + drop(beta %*% x)
  stats::dnorm(y_val, mean = mu, sd = samples$sigma, log = TRUE)
}

loglik_blasso_grid <- function(y_values, x, samples) {
  # Returns a B x K matrix whose (b,k) entry is
  # log f_{theta_b}(y_values[k] | x).
  x <- as.numeric(x)
  y_values <- as.numeric(y_values)
  beta <- as.matrix(samples$beta)
  mu <- samples$beta0 + drop(beta %*% x)
  sigma <- as.numeric(samples$sigma)
  B <- length(mu)
  K <- length(y_values)

  Z <- matrix(y_values, nrow = B, ncol = K, byrow = TRUE)
  Z <- sweep(Z, 1, mu, "-")
  Z <- sweep(Z, 1, sigma, "/")
  out <- -0.5 * Z^2
  out <- sweep(out, 1, log(sigma) + 0.5 * log(2 * pi), "-")
  out
}

loglik_train_matrix <- function(X, y, samples) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.matrix(samples$beta)
  B <- nrow(beta)
  n <- nrow(X)

  mu <- beta %*% t(X)
  mu <- sweep(mu, 1, samples$beta0, "+")
  Z <- matrix(y, nrow = B, ncol = n, byrow = TRUE) - mu
  Z <- sweep(Z, 1, samples$sigma, "/")
  out <- -0.5 * Z^2
  out <- sweep(out, 1, log(samples$sigma) + 0.5 * log(2 * pi), "-")
  out
}

# ---------------- conformal comparisons ----------------------

satisfies_quantile_condition <- function(bad_count, n, alpha) {
  # Equivalent to Fong et al.'s rank condition rank > alpha * (n + 1),
  # written in terms of bad_count = number of training predictive densities
  # that are larger than the test predictive density.
  bad_count < (1 - alpha) * (n + 1)
}

bayes_aoi_bad_count <- function(
    X,
    x_new,
    y,
    y_cand,
    samples,
    loglik_train = NULL
) {
  bayes_aoi_bad_counts_many(
    X = X,
    x_new = x_new,
    y = y,
    y_values = y_cand,
    samples = samples,
    loglik_train = loglik_train
  )[1]
}

bayes_aoi_bad_counts_many <- function(
    X,
    x_new,
    y,
    y_values,
    samples,
    loglik_train = NULL
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(loglik_train)) loglik_train <- loglik_train_matrix(X, y, samples)

  ll_new <- loglik_blasso_grid(y_values, x_new, samples) # B x K
  K <- ncol(ll_new)

  # Test posterior predictive numerator after add-one-in:
  #   sum_b f_b(y | x_new)^2.
  # The common IS denominator sum_b f_b(y | x_new) cancels in comparisons.
  log_ppd_test_num <- col_logsumexp(2 * ll_new)

  bad_count <- integer(K)
  for (k in seq_len(K)) {
    # Training predictive numerators after adding candidate test point:
    #   sum_b f_b(y | x_new) f_b(Y_i | X_i), for i=1,...,n.
    log_ppd_train_num <- col_logsumexp(sweep(loglik_train, 1, ll_new[, k], "+"))
    bad_count[k] <- sum(log_ppd_test_num[k] < log_ppd_train_num)
  }
  bad_count
}

precompute_corrected_leftout_loglik <- function(X, y, loo_samples) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  out <- vector("list", n)
  for (j in seq_len(n)) {
    out[[j]] <- loglik_blasso(y[j], X[j, ], loo_samples[[j]])
  }
  out
}

bayes_corrected_bad_count <- function(
    X,
    x_new,
    y,
    y_cand,
    loo_samples,
    leftout_loglik = NULL
) {
  bayes_corrected_bad_counts_many(
    X = X,
    x_new = x_new,
    y = y,
    y_values = y_cand,
    loo_samples = loo_samples,
    leftout_loglik = leftout_loglik
  )[1]
}

bayes_corrected_bad_counts_many <- function(
    X,
    x_new,
    y,
    y_values,
    loo_samples,
    leftout_loglik = NULL
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- length(y)
  K <- length(y_values)
  if (is.null(leftout_loglik)) {
    leftout_loglik <- precompute_corrected_leftout_loglik(X, y, loo_samples)
  }

  bad_count <- integer(K)
  for (j in seq_len(n)) {
    ll_j <- leftout_loglik[[j]]
    ll_new <- loglik_blasso_grid(y_values, x_new, loo_samples[[j]]) # B_j x K

    # Corrected/tournament comparison using proposal pi(theta | D \ {Z_j}).
    # After adding both Z_j and candidate Z_new, the common denominator is
    # sum_b f_b(Y_j | X_j) f_b(y | x_new), so it cancels. The numerators are:
    #   test:  sum_b f_b(Y_j | X_j)   f_b(y | x_new)^2
    #   train: sum_b f_b(Y_j | X_j)^2 f_b(y | x_new).
    log_ppd_test_num <- col_logsumexp(sweep(2 * ll_new, 1, ll_j, "+"))
    log_ppd_train_num <- col_logsumexp(sweep(ll_new, 1, 2 * ll_j, "+"))

    bad_count <- bad_count + as.integer(log_ppd_test_num < log_ppd_train_num)
  }
  bad_count
}

# ---------------- grid length, Fong-style ---------------------

bayes_region_grid <- function(y_grid, bad_count_function, n, alpha) {
  in_region <- logical(length(y_grid))
  for (k in seq_along(y_grid)) {
    bad <- bad_count_function(y_grid[k])
    in_region[k] <- satisfies_quantile_condition(bad, n, alpha)
  }
  in_region
}

bayes_grid_length <- function(y_grid, in_region) {
  # Fong-style grid approximation to Lebesgue length.
  if (length(y_grid) < 2) return(0)
  dy <- median(diff(y_grid))
  sum(in_region) * dy
}

bayes_grid_intervals <- function(y_grid, in_region) {
  if (!any(in_region)) return(list())
  dy <- median(diff(y_grid))
  idx <- which(in_region)
  cuts <- c(1L, which(diff(idx) > 1L) + 1L, length(idx) + 1L)
  out <- list()
  for (m in seq_len(length(cuts) - 1L)) {
    block <- idx[cuts[m]:(cuts[m + 1L] - 1L)]
    out[[length(out) + 1L]] <- c(y_grid[min(block)] - dy / 2, y_grid[max(block)] + dy / 2)
  }
  out
}

# ---------------- one split runner ----------------------------

run_bayes_lasso_conformal_split <- function(
    seed,
    X_train,
    y_train,
    X_test,
    y_test,
    alpha = 0.20,
    standardize = TRUE,
    # full posterior for Bayesian:approx
    iter = 3000,
    warmup = 1000,
    chains = 4,
    # leave-one-out posteriors for Bayesian:corrected
    loo_iter = iter,
    loo_warmup = warmup,
    loo_chains = chains,
    sigma_prior_scale = 1,
    adapt_delta = 0.8,
    max_treedepth = 10,
    # grid for length; in standardized y-scale if standardize=TRUE
    y_grid = NULL,
    grid_size = 100,
    grid_pad = 2,
    corrected = TRUE,
    max_test = NULL,
    cache_dir = NULL,
    refresh = 0,
    verbose = FALSE,
    print_progress = FALSE
) {
  check_rstan()
  if (missing(seed) || !is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("seed must be a single numeric value.", call. = FALSE)
  }
  seed <- as.integer(seed)
  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  y_train <- as.numeric(y_train)
  y_test <- as.numeric(y_test)

  if (standardize) {
    st <- standardize_split(X_train, y_train, X_test, y_test)
    X_train_s <- st$X_train
    y_train_s <- st$y_train
    X_test_s <- st$X_test
    y_test_s <- st$y_test
    y_scale <- st$y_sd
  } else {
    X_train_s <- X_train
    y_train_s <- y_train
    X_test_s <- X_test
    y_test_s <- y_test
    y_scale <- 1
  }

  n <- length(y_train_s)
  n_test <- length(y_test_s)
  if (!is.null(max_test)) n_test <- min(n_test, max_test)

  if (is.null(y_grid)) {
    lo <- min(y_train_s) - grid_pad * stats::sd(y_train_s)
    hi <- max(y_train_s) + grid_pad * stats::sd(y_train_s)
    y_grid <- seq(lo, hi, length.out = grid_size)
  }

  stan_model <- compile_fong_blasso_model()

  full_cache_file <- if (is.null(cache_dir)) NULL else file.path(cache_dir, "full_posterior.rds")
  full_samples <- read_rds_if_valid(full_cache_file)
  if (is.null(full_samples)) {
    if (verbose) cat("Fitting full posterior for Bayesian:approx\n")
    full_samples <- fit_fong_blasso_mcmc(
      X = X_train_s,
      y = y_train_s,
      stan_model = stan_model,
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed,
      sigma_prior_scale = sigma_prior_scale,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh,
      keep_fit = FALSE
    )
    save_rds_atomic(full_samples, full_cache_file)
  } else if (verbose) {
    cat("Loaded cached full posterior for Bayesian:approx\n")
  }
  full_loglik_train <- loglik_train_matrix(X_train_s, y_train_s, full_samples)

  loo_samples <- NULL
  leftout_loglik <- NULL
  if (corrected) {
    if (verbose) cat("Fitting leave-one-out posteriors for Bayesian:corrected\n")
    loo_samples <- fit_loo_fong_blasso_mcmc(
      X = X_train_s,
      y = y_train_s,
      stan_model = stan_model,
      iter = loo_iter,
      warmup = loo_warmup,
      chains = loo_chains,
      seed = seed + 10000,
      sigma_prior_scale = sigma_prior_scale,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh,
      cache_dir = if (is.null(cache_dir)) NULL else file.path(cache_dir, "loo_posteriors"),
      verbose = verbose
    )
    leftout_loglik <- precompute_corrected_leftout_loglik(X_train_s, y_train_s, loo_samples)
  }

  rows <- list()
  row_id <- 0L

  for (t in seq_len(n_test)) {
    x_new <- X_test_s[t, ]
    y_true <- y_test_s[t]
    y_values <- c(y_true, y_grid)

    approx_bad <- bayes_aoi_bad_counts_many(
      X = X_train_s,
      x_new = x_new,
      y = y_train_s,
      y_values = y_values,
      samples = full_samples,
      loglik_train = full_loglik_train
    )
    approx_bad_true <- approx_bad[1]
    approx_region <- satisfies_quantile_condition(approx_bad[-1], n, alpha)

    row_id <- row_id + 1L
    rows[[row_id]] <- data.frame(
      test_id = t,
      method = "Bayesian:approx",
      covered = satisfies_quantile_condition(approx_bad_true, n, alpha),
      length = bayes_grid_length(y_grid, approx_region),
      bad_count_true = approx_bad_true
    )

    if (corrected) {
      corrected_bad <- bayes_corrected_bad_counts_many(
        X = X_train_s,
        x_new = x_new,
        y = y_train_s,
        y_values = y_values,
        loo_samples = loo_samples,
        leftout_loglik = leftout_loglik
      )
      corrected_bad_true <- corrected_bad[1]
      corrected_region <- satisfies_quantile_condition(corrected_bad[-1], n, alpha)

      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        test_id = t,
        method = "Bayesian:corrected",
        covered = satisfies_quantile_condition(corrected_bad_true, n, alpha),
        length = bayes_grid_length(y_grid, corrected_region),
        bad_count_true = corrected_bad_true
      )
    }

    if (print_progress) {
      pct_done <- 100 * t / n_test
      cat(sprintf("Progress: %.1f%% (%d/%d test points complete)\n", pct_done, t, n_test))
      flush.console()
    }
  }

  ans <- do.call(rbind, rows)
  ans$seed <- seed
  attr(ans, "y_grid_standardized") <- y_grid
  attr(ans, "y_scale") <- y_scale
  attr(ans, "length_scale") <- "standardized_y"
  ans
}

summarize_bayes_results <- function(res) {
  aggregate(cbind(coverage = covered, length = length) ~ method, data = res, FUN = mean)
}

summarize_bayes_results_across_seeds <- function(res, alpha = 0.20) {
  if (!"seed" %in% names(res)) stop("res must contain a seed column.", call. = FALSE)
  target <- 1 - alpha
  rep_means <- stats::aggregate(
    cbind(coverage = as.numeric(covered), length = length) ~ seed + method,
    data = res,
    FUN = mean
  )

  out <- do.call(rbind, lapply(split(rep_means, rep_means$method), function(d) {
    n_reps <- nrow(d)
    cov_hat <- mean(d$coverage)
    len_hat <- mean(d$length)
    cov_se <- if (n_reps > 1L) stats::sd(d$coverage) / sqrt(n_reps) else NA_real_
    len_se <- if (n_reps > 1L) stats::sd(d$length) / sqrt(n_reps) else NA_real_
    data.frame(
      method = unique(d$method),
      coverage = cov_hat,
      coverage_se = cov_se,
      coverage_display = sprintf("%.3f (%.3f)", cov_hat, cov_se),
      coverage_not_within_3se = if (is.na(cov_se)) NA else abs(cov_hat - target) > 3 * cov_se,
      length = len_hat,
      length_se = len_se,
      length_display = sprintf("%.2f (%.2f)", len_hat, len_se),
      n_reps = n_reps
    )
  }))
  rownames(out) <- NULL
  out
}

# ---------------- diabetes helper ----------------------------

load_diabetes_lars <- function() {
  # Same 442 x 10 diabetes regression data used in the LARS/lasso literature.
  # install.packages("lars") if needed.
  if (!requireNamespace("lars", quietly = TRUE)) {
    stop("Please install lars first: install.packages('lars')")
  }
  data(diabetes, package = "lars")
  list(X = as.matrix(diabetes$x), y = as.numeric(diabetes$y))
}

make_train_test_split <- function(X, y, train_frac = 0.7, seed = 100) {
  set.seed(seed)
  n <- length(y)
  n_train <- floor(train_frac * n)
  idx <- s ample.int(n)
  train_idx <- idx[seq_len(n_train)]
  test_idx <- idx[(n_train + 1):n]
  list(
    X_train = X[train_idx, , drop = FALSE],
    y_train = y[train_idx],
    X_test = X[test_idx, , drop = FALSE],
    y_test = y[test_idx],
    train_idx = train_idx,
    test_idx = test_idx
  )
}
