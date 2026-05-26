# ============================================================
# Stability frontiers for the four approximation families in the paper
#
# Notation:
#   - s(z ; D) is the score obtained by training on the dataset D
#     and evaluating the score at z.
#   - s_approx is the approximate score from the paper.
#
# For each approximation family, we plot the empirical epsilon(nu)
# frontier, where epsilon(nu) is the empirical (1 - nu)-quantile of
# the corresponding score discrepancy Delta.
#
# Approximate-method discrepancy:
#   Delta_i^{approx}
#   = | s_approx(Z_i ; D_{n+1}\setminus{Z_i} + {Z_i})
#       - s(Z_i ; D_{n+1}) |.
#
# Corrected-method discrepancy:
#   Delta_i^{corr}
#   = | s_approx(Z_i ; D_{n\setminus i} + {Z_i, Z_{n+1}})
#       - s_approx(Z_i ; D_{n+1}\setminus{Z_i} + {Z_i}) |.
#
# By exchangeability, these are the score discrepancies corresponding
# to the stability comparisons used in the paper.
#
# Approximation families:
#   1) Deletion
#   2) Rounding
#   3) One-step update
#   4) Bayesian AOI
# ============================================================

source("~/Generalization of cross-conformal/conformal_simulation_files/methods.R")

# ---------- Small helpers ----------
abs_residual_score <- function(x, y, beta) {
  abs(y - ols_predict(matrix(x, nrow = 1), beta))
}

nearest_grid_value <- function(y, grid) {
  grid[which.min(abs(grid - y))]
}

# Smallest empirical epsilon such that P(Delta <= epsilon) >= level.
quantile_higher <- function(x, level) {
  x <- sort(as.numeric(x))
  n <- length(x)

  if (level <= 0) return(x[1])
  if (level >= 1) return(x[n])

  idx <- ceiling(level * n)
  idx <- max(1L, min(n, idx))
  x[idx]
}

# Monte Carlo approximation of E[exp(log_w)] computed stably from log_w.
mc_meanexp <- function(log_w) {
  exp(logsumexp(log_w) - log(length(log_w)))
}

# ---------- One Monte Carlo replicate ----------
simulate_one_replicate_all4 <- function(
    n,
    p,
    sigma,
    eta,
    M,
    T = 100,
    mu_prior = 10,
    sig_prior = 1,
    sig_lik_model = 1
) {
  beta_true <- rand_beta_unit_scaled(p)

  X_all <- matrix(rnorm((n + 1) * p), nrow = n + 1, ncol = p)
  y_all <- drop(X_all %*% beta_true + sigma * rnorm(n + 1))

  X_train <- X_all[1:n, , drop = FALSE]
  y_train <- y_all[1:n]
  x_new <- X_all[n + 1, ]
  y_new <- y_all[n + 1]

  grid <- make_grid(y_train, M)
  ybar_train <- vapply(y_train, nearest_grid_value, numeric(1), grid = grid)
  ybar_new <- nearest_grid_value(y_new, grid)

  # Reference score s(Z_i ; D_{n+1}) for the first three families.
  beta_full <- ols_beta(X_all, y_all)
  score_full_all <- abs(y_train - ols_predict(X_train, beta_full))

  # Posterior draws needed for the Bayesian family.
  mu_0 <- rep(mu_prior, p)
  Sigma_0 <- diag(sig_prior^2, p)

  theta_full <- bayesian_lr(
    mu_0 = mu_0,
    Sigma_0 = Sigma_0,
    X = X_all,
    y = y_all,
    T = T,
    sig_lik_model = sig_lik_model
  )

  theta_full_loo <- theta_loo(
    mu_0 = mu_0,
    Sigma_0 = Sigma_0,
    X = X_all,
    y = y_all,
    T = T,
    sig_lik_model = sig_lik_model
  )

  theta_train_loo <- theta_loo(
    mu_0 = mu_0,
    Sigma_0 = Sigma_0,
    X = X_train,
    y = y_train,
    T = T,
    sig_lik_model = sig_lik_model
  )

  delta_deletion_approx <- numeric(n)
  delta_rounding_approx <- numeric(n)
  delta_onestep_approx <- numeric(n)
  delta_bayesian_approx <- numeric(n)

  delta_deletion_corrected <- numeric(n)
  delta_rounding_corrected <- numeric(n)
  delta_onestep_corrected <- numeric(n)
  delta_bayesian_corrected <- numeric(n)

  for (i in seq_len(n)) {
    xi <- X_train[i, ]
    yi <- y_train[i]
    yi_bar <- ybar_train[i]

    X_minus_i <- X_train[-i, , drop = FALSE]
    y_minus_i <- y_train[-i]

    # D_{n+1} \ {Z_i} = (D_n \ {Z_i}) U {Z_{n+1}}
    X_full_minus_i <- rbind(X_minus_i, matrix(x_new, nrow = 1))
    y_full_minus_i <- c(y_minus_i, y_new)

    score_full_i <- score_full_all[i]

    # --------------------------------------------------
    # Deletion family
    # --------------------------------------------------
    # s_approx(Z_i ; D_{n+1}\setminus{Z_i} + {Z_i})
    # = s(Z_i ; D_{n+1}\setminus{Z_i})
    beta_one_special_delete <- ols_beta(X_full_minus_i, y_full_minus_i)
    score_one_special_delete <- abs_residual_score(xi, yi, beta_one_special_delete)

    # s_approx(Z_i ; D_{n\setminus i} + {Z_i, Z_{n+1}})
    # = s(Z_i ; D_{n\setminus i})
    beta_two_special_delete <- ols_beta(X_minus_i, y_minus_i)
    score_two_special_delete <- abs_residual_score(xi, yi, beta_two_special_delete)

    delta_deletion_approx[i] <- abs(score_one_special_delete - score_full_i)
    delta_deletion_corrected[i] <- abs(score_two_special_delete - score_one_special_delete)

    # --------------------------------------------------
    # Rounding family
    # --------------------------------------------------
    # s_approx(Z_i ; D_{n+1}\setminus{Z_i} + {Z_i})
    # = score at Z_i when only Z_i is rounded.
    y_one_special_round <- y_all
    y_one_special_round[i] <- yi_bar
    beta_one_special_round <- ols_beta(X_all, y_one_special_round)
    score_one_special_round <- abs_residual_score(xi, yi, beta_one_special_round)

    # s_approx(Z_i ; D_{n\setminus i} + {Z_i, Z_{n+1}})
    # = score at Z_i when both Z_i and Z_{n+1} are rounded.
    X_two_special_round <- rbind(
      X_minus_i,
      matrix(xi, nrow = 1),
      matrix(x_new, nrow = 1)
    )
    y_two_special_round <- c(y_minus_i, yi_bar, ybar_new)
    beta_two_special_round <- ols_beta(X_two_special_round, y_two_special_round)
    score_two_special_round <- abs_residual_score(xi, yi, beta_two_special_round)

    delta_rounding_approx[i] <- abs(score_one_special_round - score_full_i)
    delta_rounding_corrected[i] <- abs(score_two_special_round - score_one_special_round)

    # --------------------------------------------------
    # One-step update family
    # --------------------------------------------------
    # s_approx(Z_i ; D_{n+1}\setminus{Z_i} + {Z_i})
    # = score at Z_i after a one-step update that adds back Z_i.
    g_i <- eta / (n + 1) * sum(xi * xi)
    resid_one_special <- yi - ols_predict(matrix(xi, nrow = 1), beta_one_special_delete)
    score_one_special_onestep <- abs((1 - g_i) * resid_one_special)

    # s_approx(Z_i ; D_{n\setminus i} + {Z_i, Z_{n+1}})
    # = score at Z_i after a one-step update using both Z_i and Z_{n+1}.
    a_i <- eta / (n + 1) * sum(xi * x_new)
    resid_two_special_i <- yi - ols_predict(matrix(xi, nrow = 1), beta_two_special_delete)
    resid_two_special_new <- y_new - ols_predict(matrix(x_new, nrow = 1), beta_two_special_delete)
    score_two_special_onestep <- abs((1 - g_i) * resid_two_special_i - a_i * resid_two_special_new)

    delta_onestep_approx[i] <- abs(score_one_special_onestep - score_full_i)
    delta_onestep_corrected[i] <- abs(score_two_special_onestep - score_one_special_onestep)

    # --------------------------------------------------
    # Bayesian AOI family
    # --------------------------------------------------
    # Here s(z ; D) is the full posterior predictive score and
    # s_approx is the AOI / corrected Bayesian score from the paper.

    theta_full_minus_i <- theta_full_loo[i, , ]
    theta_train_minus_i <- theta_train_loo[i, , ]

    loglik_full_i <- loglik_bayes(
      yi, xi, theta_full,
      sig_lik_model = sig_lik_model
    )
    loglik_full_minus_i_i <- loglik_bayes(
      yi, xi, theta_full_minus_i,
      sig_lik_model = sig_lik_model
    )
    loglik_train_minus_i_i <- loglik_bayes(
      yi, xi, theta_train_minus_i,
      sig_lik_model = sig_lik_model
    )
    loglik_train_minus_i_new <- loglik_bayes(
      y_new, x_new, theta_train_minus_i,
      sig_lik_model = sig_lik_model
    )

    # s(Z_i ; D_{n+1})
    score_full_bayes <- mc_meanexp(loglik_full_i)

    # s_approx(Z_i ; D_{n+1}\setminus{Z_i} + {Z_i})
    score_one_special_bayes <- mc_meanexp(2.0 * loglik_full_minus_i_i)

    # s_approx(Z_i ; D_{n\setminus i} + {Z_i, Z_{n+1}})
    score_two_special_bayes <- mc_meanexp(2.0 * loglik_train_minus_i_i + loglik_train_minus_i_new)

    delta_bayesian_approx[i] <- abs(score_one_special_bayes - score_full_bayes)
    delta_bayesian_corrected[i] <- abs(score_two_special_bayes - score_one_special_bayes)
  }

  data.frame(
    family = rep(c(
      "Example 0: Deletion",
      "Example 1: Rounding",
      "Example 2: One-step update",
      "Example 3: Bayesian"
    ), each = 2 * n),
    version = rep(rep(c("Approximate", "Corrected"), each = n), times = 4),
    delta = c(
      delta_deletion_approx,
      delta_deletion_corrected,
      delta_rounding_approx,
      delta_rounding_corrected,
      delta_onestep_approx,
      delta_onestep_corrected,
      delta_bayesian_approx,
      delta_bayesian_corrected
    )
  )
}

# ---------- Run many replicates ----------
simulate_stability_frontiers_all4 <- function(
    n = 100,
    p = 20,
    reps = 500,
    sigma = 1,
    eta = 0.1,
    M = 10,
    T = 100,
    mu_prior = 10,
    sig_prior = 1,
    sig_lik_model = 1,
    seed = 123,
    show_progress = TRUE
) {
  set.seed(seed)

  out <- vector("list", reps)

  for (rep in seq_len(reps)) {
    out[[rep]] <- simulate_one_replicate_all4(
      n = n,
      p = p,
      sigma = sigma,
      eta = eta,
      M = M,
      T = T,
      mu_prior = mu_prior,
      sig_prior = sig_prior,
      sig_lik_model = sig_lik_model
    )

    if (show_progress && (rep %% 10 == 0 || rep == reps)) {
      message(sprintf("Completed %d / %d replicates", rep, reps))
    }
  }

  do.call(rbind, out)
}

# ---------- Convert deltas into epsilon(nu) frontiers ----------
build_frontier_df_all4 <- function(
    delta_df,
    nu_grid = seq(0.001, 0.999, length.out = 300)
) {
  groups <- unique(delta_df[, c("family", "version")])
  out <- vector("list", nrow(groups))

  for (j in seq_len(nrow(groups))) {
    fam <- groups$family[j]
    ver <- groups$version[j]

    vals <- delta_df$delta[
      delta_df$family == fam & delta_df$version == ver
    ]

    out[[j]] <- data.frame(
      nu = nu_grid,
      epsilon = vapply(1 - nu_grid, function(level) {
        quantile_higher(vals, level)
      }, numeric(1)),
      family = fam,
      version = ver
    )
  }

  frontier_df <- do.call(rbind, out)

  frontier_df$family <- factor(
    frontier_df$family,
    levels = c(
      "Example 0: Deletion",
      "Example 1: Rounding",
      "Example 2: One-step update",
      "Example 3: Bayesian"
    )
  )

  frontier_df$version <- factor(
    frontier_df$version,
    levels = c("Approximate", "Corrected")
  )

  frontier_df
}

# ---------- Plot ----------
plot_frontier_all4 <- function(frontier_df, n, p, eta, M, T) {
  method_cols <- c(
    "Approximate" = "#2C7FB8",
    "Corrected"   = "#D7191C"
  )
  
  ymax <- max(frontier_df$epsilon, na.rm = TRUE)
  
  # One fake x-axis title per facet
  xlab_df <- data.frame(
    family = factor(levels(frontier_df$family), levels = levels(frontier_df$family)),
    nu = 0.5,
    epsilon = -0.105 * ymax,
    label = "nu"
  )
  
  ggplot2::ggplot(
    frontier_df,
    ggplot2::aes(
      x = nu,
      y = epsilon,
      color = version,
      group = interaction(family, version)
    )
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_text(
      data = xlab_df,
      ggplot2::aes(x = nu, y = epsilon, label = label),
      inherit.aes = FALSE,
      parse = TRUE,
      size = 12
    ) +
    ggplot2::facet_wrap(~ family, nrow = 1) +
    ggplot2::scale_color_manual(
      name = "Method",
      values = method_cols,
      breaks = c("Approximate", "Corrected"),
      labels = c("Approximate", "Tournament"),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, ymax),
      clip = "off"
    ) +
    ggplot2::labs(
      x = NULL,
      y = expression(epsilon)
    ) +
    ggplot2::theme_bw(base_size = 26) +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      strip.background = ggplot2::element_rect(fill = "grey85"),
      strip.text = ggplot2::element_text(face = "bold", size = 24),
      axis.title.x = ggplot2::element_text(size = 36),
      axis.title.y = ggplot2::element_text(size = 36),
      axis.text = ggplot2::element_text(size = 26),
      legend.title = ggplot2::element_text(size = 26),
      legend.text = ggplot2::element_text(size = 26),
      
      # Extra bottom space for the four ν labels
      plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 40, l = 5.5)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(linewidth = 1.4)
      )
    )
}
# ---------- Convenience wrapper ----------
run_stability_plot_all4 <- function(
    n = 100,
    p = 20,
    reps = 200,
    sigma = 1,
    eta = 0.1,
    M = 10,
    T = 100,
    mu_prior = 10,
    sig_prior = 1,
    sig_lik_model = 1,
    seed = 123,
    show_progress = TRUE
) {
  all_deltas <- simulate_stability_frontiers_all4(
    n = n,
    p = p,
    reps = reps,
    sigma = sigma,
    eta = eta,
    M = M,
    T = T,
    mu_prior = mu_prior,
    sig_prior = sig_prior,
    sig_lik_model = sig_lik_model,
    seed = seed,
    show_progress = show_progress
  )

  frontier_df <- build_frontier_df_all4(all_deltas)

  p_frontier <- plot_frontier_all4(
    frontier_df = frontier_df,
    n = n,
    p = p,
    eta = eta,
    M = M,
    T = T
  )
  list(
    deltas = all_deltas,
    frontier = frontier_df,
    plot = p_frontier
  )
}



res_stable <- run_stability_plot_all4(
  n = 100,
  p = 80,
  reps = 200,
  sigma = 1,
  eta = 0.1,
  M = 10,
  T = 100,
  seed = 123
)
print(res_stable$plot)

