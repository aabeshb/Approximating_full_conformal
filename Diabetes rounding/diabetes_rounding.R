# ============================================================
# Rounding implementation for the Diabetes data:
#   1) Full conformal prediction with the same grid-search idea
#      used in Lee-Zhang's real-data code, using Huber-RLM fits.
#   2) Tournament-corrected rounding-one-at-a-time.
#
# ===========================================================

# ---------- basic helpers ----------
pop_sd <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x)
  sqrt(mean((x - mu)^2))
}

np_quantile_higher <- function(x, prob) {
  x <- sort(as.numeric(x))
  n <- length(x)
  if (n == 0L) stop("x must be non-empty")
  prob <- min(max(prob, 0), 1)
  idx <- 1L + ceiling((n - 1L) * prob)
  idx <- max(1L, min(n, idx))
  x[idx]
}

nearest_grid_index <- function(v, grid) {
  which.min(abs(grid - v))
}

merge_intervals <- function(intervals) {
  if (length(intervals) == 0L) return(list())

  M <- do.call(rbind, intervals)
  M <- M[order(M[, 1], M[, 2]), , drop = FALSE]

  out <- list()
  curL <- M[1, 1]
  curR <- M[1, 2]

  if (nrow(M) >= 2L) {
    for (i in 2:nrow(M)) {
      L <- M[i, 1]
      R <- M[i, 2]
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

interval_union_length <- function(intervals) {
  ints <- merge_intervals(intervals)
  if (length(ints) == 0L) return(0)

  total <- 0
  for (int in ints) {
    total <- total + (int[2] - int[1])
  }
  total
}

interval_contains <- function(intervals, y) {
  if (length(intervals) == 0L) return(FALSE)
  for (int in intervals) {
    if (y >= int[1] && y <= int[2]) return(TRUE)
  }
  FALSE
}

# ---------- Diabetes data ----------
# Uses the diabetes data from the lars package.
# We then apply the same preprocessing convention used in the paper's code:
#   X <- (X - mean(X_j)) / pop_sd(X_j) / sqrt(p)
#   y <- (y - mean(y)) / pop_sd(y)
load_diabetes_data <- function() {
  if (!requireNamespace("lars", quietly = TRUE)) {
    stop("Please install the lars package: install.packages('lars')")
  }

  env <- new.env(parent = emptyenv())
  utils::data("diabetes", package = "lars", envir = env)
  ds <- env$diabetes

  X <- as.matrix(ds$x)
  y <- as.numeric(ds$y)

  X <- sweep(X, 2, colMeans(X), "-")
  sx <- apply(X, 2, pop_sd)
  sx[sx == 0] <- 1
  X <- sweep(X, 2, sx, "/")
  X <- X / sqrt(ncol(X))

  y <- y - mean(y)
  sy <- pop_sd(y)
  if (sy == 0) sy <- 1
  y <- y / sy

  list(X = X, y = y)
}

# ---------- train/test split ----------
train_test_split_r <- function(X, y, test_size) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)

  if (length(test_size) != 1L || test_size < 1L || test_size >= n) {
    stop("test_size must be an integer in {1, ..., n-1}")
  }

  idx_test <- sample.int(n, size = as.integer(test_size), replace = FALSE)
  idx_train <- setdiff(seq_len(n), idx_test)

  list(
    train = list(X = X[idx_train, , drop = FALSE], y = y[idx_train]),
    test  = list(X = X[idx_test,  , drop = FALSE], y = y[idx_test]),
    idx_train = idx_train,
    idx_test = idx_test
  )
}

# ---------- common candidate grid ----------
# Matches the paper's grid-search style:
#   seq(min(y_train) - sd(y_train), max(y_train) + sd(y_train), length.out = 100)
# using population sd.
make_repo_grid <- function(y, grid_size = 100L) {
  y <- as.numeric(y)
  seq(min(y) - pop_sd(y), max(y) + pop_sd(y), length.out = as.integer(grid_size))
}

# ---------- Huber-RLM fit ----------
# Objective:
#   mean(huber(y - X beta, thres))/2 + lamb/2 * ||beta||_2^2
# with no intercept.
huber_value <- function(r, thres) {
  ifelse(abs(r) <= thres, r^2, 2 * thres * abs(r) - thres^2)
}

huber_grad <- function(r, thres) {
  ifelse(abs(r) <= thres, r, thres * sign(r))
}

fit_rlm_huber <- function(X, y,
                          thres = 1,
                          lamb = 2,
                          start = NULL,
                          maxit = 1000,
                          reltol = 1e-9) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(start)) {
    start <- rep(0, p)
  } else {
    start <- as.numeric(start)
    if (length(start) != p || anyNA(start)) {
      start <- rep(0, p)
    }
  }

  obj_fn <- function(beta) {
    r <- as.numeric(y - X %*% beta)
    mean(huber_value(r, thres)) / 2 + (lamb / 2) * sum(beta^2)
  }

  grad_fn <- function(beta) {
    r <- as.numeric(y - X %*% beta)
    psi <- huber_grad(r, thres)
    as.numeric(-(crossprod(X, psi) / n) + lamb * beta)
  }

  fit <- stats::optim(
    par = start,
    fn = obj_fn,
    gr = grad_fn,
    method = "BFGS",
    control = list(reltol = reltol, maxit = maxit)
  )

  beta_hat <- as.numeric(fit$par)
  if (length(beta_hat) != p || anyNA(beta_hat) || any(!is.finite(beta_hat))) {
    stop("optim() failed to return a valid coefficient vector.")
  }
  beta_hat
}

predict_rlm_huber <- function(X, beta) {
  as.numeric(as.matrix(X) %*% as.numeric(beta))
}

# ============================================================
# Full conformal prediction using the gridding approach
# ============================================================
# This is the exact score s(z; D): fit the Huber-RLM on the augmented data
# and evaluate the absolute residual score at every point.

fullcp_interval_rlm <- function(X, y, x_new,
                                alpha = 0.1,
                                thres = 1,
                                lamb = 2,
                                grid = NULL,
                                grid_size = 100L,
                                maxit = 1000,
                                reltol = 1e-9) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  x_new <- as.numeric(x_new)

  if (is.null(grid)) grid <- make_repo_grid(y, grid_size = grid_size)
  M <- length(grid)
  X_aug <- rbind(X, matrix(x_new, nrow = 1L))
  p <- ncol(X)

  lower <- NA_real_
  upper <- NA_real_

  beta_start <- rep(0, p)
  for (m in seq_len(M)) {
    z <- grid[m]
    y_aug <- c(y, z)
    beta_hat <- fit_rlm_huber(
      X = X_aug,
      y = y_aug,
      thres = thres,
      lamb = lamb,
      start = beta_start,
      maxit = maxit,
      reltol = reltol
    )
    beta_start <- beta_hat

    pred <- predict_rlm_huber(X_aug, beta_hat)
    S <- abs(pred - y_aug)
    Q <- np_quantile_higher(S, 1 - alpha)

    if (S[length(S)] <= Q) {
      lower <- z
      break
    }
  }

  if (!is.finite(lower)) {
    return(c(NA_real_, NA_real_))
  }

  beta_start <- rep(0, p)
  for (m in M:1L) {
    z <- grid[m]
    y_aug <- c(y, z)
    beta_hat <- fit_rlm_huber(
      X = X_aug,
      y = y_aug,
      thres = thres,
      lamb = lamb,
      start = beta_start,
      maxit = maxit,
      reltol = reltol
    )
    beta_start <- beta_hat

    pred <- predict_rlm_huber(X_aug, beta_hat)
    S <- abs(pred - y_aug)
    Q <- np_quantile_higher(S, 1 - alpha)

    if (S[length(S)] <= Q) {
      upper <- z
      break
    }
  }

  c(lower, upper)
}

# ============================================================
# Tournament-corrected rounding-one-at-a-time
# ============================================================
# This is the corrected score s_approx(z; D_{n\setminus i} + {z, z'}).
# For each i, only the i-th training response is rounded to the common grid.
# The score remains the absolute residual evaluated at the point of interest.

round_one_intervals_rlm <- function(X, y, x_new,
                                    alpha = 0.1,
                                    thres = 1,
                                    lamb = 2,
                                    grid = NULL,
                                    grid_size = 100L,
                                    maxit = 1000,
                                    reltol = 1e-9) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  x_new <- as.numeric(x_new)
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(grid)) grid <- make_repo_grid(y, grid_size = grid_size)
  rounded_y <- grid[vapply(y, nearest_grid_index, integer(1L), grid = grid)]

  mids <- (grid[-length(grid)] + grid[-1L]) / 2
  left_edges  <- c(-Inf, mids)
  right_edges <- c(mids, Inf)

  tau <- (1 - alpha) * (n + 1)
  K <- n - floor(tau)

  X_aug <- rbind(X, matrix(x_new, nrow = 1L))
  intervals <- list()

  for (m in seq_along(grid)) {
    y_m <- grid[m]
    evx <- numeric(0)
    evd <- integer(0)
    beta_start <- rep(0, p)

    for (i in seq_len(n)) {
      y_mod <- y
      y_mod[i] <- rounded_y[i]
      y_aug <- c(y_mod, y_m)

      beta_hat <- fit_rlm_huber(
        X = X_aug,
        y = y_aug,
        thres = thres,
        lamb = lamb,
        start = beta_start,
        maxit = maxit,
        reltol = reltol
      )
      beta_start <- beta_hat

      pred_all <- predict_rlm_huber(X_aug, beta_hat)
      mu_i <- pred_all[n + 1L]
      r_i <- abs(y[i] - pred_all[i])

      L <- mu_i - r_i
      R <- mu_i + r_i
      if (is.finite(L) && is.finite(R) && L < R) {
        evx <- c(evx, L, R)
        evd <- c(evd, +1L, -1L)
      }
    }

    if (length(evx) == 0L) next

    events_df <- stats::aggregate(delta ~ x, data = data.frame(x = evx, delta = evd), FUN = sum)
    events_df <- events_df[order(events_df$x), , drop = FALSE]

    coverage <- 0L
    prev_x <- NULL

    for (k in seq_len(nrow(events_df))) {
      xk <- events_df$x[k]
      if (!is.null(prev_x) && coverage >= K && xk > prev_x) {
        Lc <- max(prev_x, left_edges[m])
        Rc <- min(xk, right_edges[m])
        if (Lc < Rc) {
          intervals[[length(intervals) + 1L]] <- c(Lc, Rc)
        }
      }
      coverage <- coverage + events_df$delta[k]
      prev_x <- xk
    }
  }

  merge_intervals(intervals)
}

# ---------- one-point evaluation ----------
evaluate_one_test_point <- function(X, y, x_new, y_true,
                                    alpha = 0.1,
                                    thres = 1,
                                    lamb = 2,
                                    grid = NULL,
                                    grid_size = 100L,
                                    maxit = 1000,
                                    reltol = 1e-9) {
  if (is.null(grid)) grid <- make_repo_grid(y, grid_size = grid_size)

  full_int <- fullcp_interval_rlm(
    X = X,
    y = y,
    x_new = x_new,
    alpha = alpha,
    thres = thres,
    lamb = lamb,
    grid = grid,
    maxit = maxit,
    reltol = reltol
  )

  corrected_ints <- round_one_intervals_rlm(
    X = X,
    y = y,
    x_new = x_new,
    alpha = alpha,
    thres = thres,
    lamb = lamb,
    grid = grid,
    maxit = maxit,
    reltol = reltol
  )

  full_cover <- if (anyNA(full_int)) {
    NA
  } else {
    (y_true >= full_int[1]) && (y_true <= full_int[2])
  }

  full_length <- if (anyNA(full_int)) {
    NA_real_
  } else {
    max(0, full_int[2] - full_int[1])
  }

  corrected_cover <- interval_contains(corrected_ints, y_true)
  corrected_length <- interval_union_length(corrected_ints)

  list(
    full_interval = full_int,
    corrected_intervals = corrected_ints,
    summary = data.frame(
      y_true = y_true,
      full_cover = full_cover,
      full_length = full_length,
      corrected_cover = corrected_cover,
      corrected_length = corrected_length
    )
  )
}

# ---------- repeated experiment ----------
run_diabetes_experiment <- function(n_iter = 100L,
                                    test_size = 1L,
                                    alpha = 0.1,
                                    seed = 1L,
                                    thres = 1,
                                    lamb = 2,
                                    grid_size = 100L,
                                    maxit = 1000,
                                    reltol = 1e-9,
                                    verbose = TRUE,
                                    print_progress = FALSE,
                                    progress_every = 1L) {
  dat <- load_diabetes_data()
  set.seed(seed)

  progress_every <- as.integer(progress_every)
  if (!is.finite(progress_every) || is.na(progress_every) || progress_every < 1L) {
    progress_every <- 1L
  }

  iteration_results <- data.frame(
    iter = seq_len(n_iter),
    n_test = NA_integer_,
    full_coverage = NA_real_,
    full_avg_length = NA_real_,
    corrected_coverage = NA_real_,
    corrected_avg_length = NA_real_
  )

  point_results <- vector("list", n_iter)

  for (b in seq_len(n_iter)) {
    if (isTRUE(verbose) && ((b <= 5L) || (b %% 10L == 0L))) {
      message("iteration ", b, " / ", n_iter)
    }

    split <- train_test_split_r(dat$X, dat$y, test_size = test_size)
    X_train <- split$train$X
    y_train <- split$train$y
    X_test <- split$test$X
    y_test <- split$test$y
    grid <- make_repo_grid(y_train, grid_size = grid_size)

    m <- nrow(X_test)
    point_df <- data.frame(
      iter = rep(b, m),
      test_id = seq_len(m),
      y_true = y_test,
      full_cover = NA,
      full_length = NA_real_,
      corrected_cover = NA,
      corrected_length = NA_real_
    )

    for (j in seq_len(m)) {
      out_j <- evaluate_one_test_point(
        X = X_train,
        y = y_train,
        x_new = X_test[j, ],
        y_true = y_test[j],
        alpha = alpha,
        thres = thres,
        lamb = lamb,
        grid = grid,
        maxit = maxit,
        reltol = reltol
      )

      point_df$full_cover[j] <- out_j$summary$full_cover
      point_df$full_length[j] <- out_j$summary$full_length
      point_df$corrected_cover[j] <- out_j$summary$corrected_cover
      point_df$corrected_length[j] <- out_j$summary$corrected_length

      if (isTRUE(print_progress) && ((j %% progress_every == 0L) || (j == 1L) || (j == m))) {
        full_so_far <- if (all(is.na(point_df$full_cover[seq_len(j)]))) {
          NA_real_
        } else {
          100 * mean(point_df$full_cover[seq_len(j)], na.rm = TRUE)
        }
        corrected_so_far <- if (all(is.na(point_df$corrected_cover[seq_len(j)]))) {
          NA_real_
        } else {
          100 * mean(point_df$corrected_cover[seq_len(j)], na.rm = TRUE)
        }
        pct_done <- 100 * j / m
        message(sprintf(
          "iteration %d / %d | test point %d / %d (%.1f%% complete) | FullCP covered so far: %s | Corrected covered so far: %s",
          b,
          n_iter,
          j,
          m,
          pct_done,
          if (is.na(full_so_far)) "NA" else sprintf("%.1f%%", full_so_far),
          if (is.na(corrected_so_far)) "NA" else sprintf("%.1f%%", corrected_so_far)
        ))
      }
    }

    point_results[[b]] <- point_df

    iteration_results$n_test[b] <- m
    iteration_results$full_coverage[b] <- mean(point_df$full_cover, na.rm = TRUE)
    iteration_results$full_avg_length[b] <- mean(point_df$full_length, na.rm = TRUE)
    iteration_results$corrected_coverage[b] <- mean(point_df$corrected_cover, na.rm = TRUE)
    iteration_results$corrected_avg_length[b] <- mean(point_df$corrected_length, na.rm = TRUE)
  }

  point_results <- do.call(rbind, point_results)

  summary <- data.frame(
    method = c("FullCP_RLM_grid", "RoundOne_corrected_RLM"),
    mean_coverage = c(
      mean(iteration_results$full_coverage, na.rm = TRUE),
      mean(iteration_results$corrected_coverage, na.rm = TRUE)
    ),
    se_coverage = c(
      stats::sd(iteration_results$full_coverage, na.rm = TRUE) / sqrt(n_iter),
      stats::sd(iteration_results$corrected_coverage, na.rm = TRUE) / sqrt(n_iter)
    ),
    mean_length = c(
      mean(iteration_results$full_avg_length, na.rm = TRUE),
      mean(iteration_results$corrected_avg_length, na.rm = TRUE)
    ),
    se_length = c(
      stats::sd(iteration_results$full_avg_length, na.rm = TRUE) / sqrt(n_iter),
      stats::sd(iteration_results$corrected_avg_length, na.rm = TRUE) / sqrt(n_iter)
    )
  )

  list(
    iteration_results = iteration_results,
    point_results = point_results,
    summary = summary,
    settings = list(
      n_iter = n_iter,
      test_size = test_size,
      alpha = alpha,
      seed = seed,
      thres = thres,
      lamb = lamb,
      grid_size = grid_size,
      maxit = maxit,
      reltol = reltol,
      verbose = verbose,
      print_progress = print_progress,
      progress_every = progress_every
    )
  )
}
