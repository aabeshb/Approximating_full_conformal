# ============================================================
# Simulation and plots for approximate and corrected conformal methods
# ============================================================
# This script uses the same linear-model data-generating process for all
# four method families:
#   1. Deletion
#   2. Rounding
#   3. One-step update
#   4. Bayesian AOI-PPD
#
# Put this file and methods.R in the same directory, then run this script.

## This will produce the coverage and length plot in Section 4.1.
# ============================================================

library(dplyr)
library(ggplot2)

# Source methods.R from the same directory as this script whenever possible.
get_script_dir <- function() {
  # Works when called as: Rscript /path/to/plot.R
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_idx <- grep(paste0("^", file_arg), cmd_args)
  if (length(file_idx) > 0L) {
    script_path <- sub(paste0("^", file_arg), "", cmd_args[file_idx[1L]])
    if (nzchar(script_path)) {
      return(dirname(normalizePath(script_path)))
    }
  }

  # Works in many source() calls.
  frame_file <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(frame_file) && nzchar(frame_file)) {
    return(dirname(normalizePath(frame_file)))
  }

  # Works when the script is opened and run from RStudio.
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    rstudio_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                             error = function(e) "")
    if (nzchar(rstudio_path)) {
      return(dirname(normalizePath(rstudio_path)))
    }
  }

  getwd()
}

script_dir <- get_script_dir()
source(file.path(script_dir, "methods.R"))

# ============================================================
# Simulation
# ============================================================

simulate_dim_curve_all <- function(
    n = 100,
    p_values = seq(20, 100, 5),
    reps = 100,
    alpha = 0.10,
    M = 10,
    sigma = 1,
    seed = 1,
    eta_onegrad = 10,
    eta_corrected_onegrad = 10,
    T_bayes = 100,
    mu_prior = 10.0,
    sig_prior = 1.0,
    sig_lik_model = 1.0,
    bayes_min_y_aoi = -80,
    bayes_max_y_aoi = 80,
    bayes_min_y_corrected = -80,
    bayes_max_y_corrected = 80,
    bayes_length_grid_size = 0.1,
    bayes_shrink = 10,
    include_bayesian = TRUE
) {
  set.seed(seed)

  total_iters <- length(p_values) * reps
  iter <- 0L

  rows_per_iter <- if (include_bayesian) 8L else 6L
  out <- vector("list", total_iters * rows_per_iter)
  idx <- 1L

  clean_length <- function(x) {
    if (!is.finite(x)) return(NA_real_)
    x
  }

  for (p in p_values) {
    for (b in seq_len(reps)) {
      iter <- iter + 1L
      if (iter %% max(1L, floor(total_iters / 100)) == 0L || iter == total_iters) {
        pct <- 100 * iter / total_iters
        cat(sprintf("\rProgress: %d/%d (%.1f%%)", iter, total_iters, pct))
        flush.console()
      }

      # Shared data generation for all methods.
      X <- matrix(rnorm(n * p), nrow = n, ncol = p)
      x_new <- rnorm(p)
      beta_true <- rand_beta_unit_scaled(p, scale = sqrt(10))

      y <- drop(X %*% beta_true + rnorm(n, sd = sigma))
      y_new <- drop(sum(x_new * beta_true) + rnorm(1, sd = sigma))

      grid_vals <- make_grid(y, M)

      # ------------------------------------------------------
      # Example 0: Deletion
      # ------------------------------------------------------
      cov_deletion_approx <- naive_cp_contains_linear_abs(X, y, x_new, y_new, alpha)
      len_deletion_approx <- naive_cp_length_linear_abs(X, y, x_new, alpha)

      cov_deletion_corrected <- loo_cross_conformal_contains_linear_abs(X, y, x_new, y_new, alpha)
      len_deletion_corrected <- loo_cross_conformal_length_linear_abs(X, y, x_new, alpha)

      out[[idx]] <- data.frame(
        rep = b,
        p = p,
        method = "Deletion:approx",
        covered = cov_deletion_approx,
        length = clean_length(len_deletion_approx)
      )
      idx <- idx + 1L

      out[[idx]] <- data.frame(
        rep = b,
        p = p,
        method = "Deletion:corrected",
        covered = cov_deletion_corrected,
        length = clean_length(len_deletion_corrected)
      )
      idx <- idx + 1L

      # ------------------------------------------------------
      # Example 1: Rounding
      # ------------------------------------------------------
      cov_rounding_approx <- gridcp_contains(X, y, x_new, y_new, grid_vals, alpha)
      len_rounding_approx <- gridcp_length(X, y, x_new, grid_vals, alpha)

      cov_rounding_corrected <- round_one_contains(X, y, x_new, y_new, grid_vals, alpha)
      len_rounding_corrected <- round_one_length(X, y, x_new, grid_vals, alpha)

      out[[idx]] <- data.frame(
        rep = b,
        p = p,
        method = "Rounding:approx",
        covered = cov_rounding_approx,
        length = clean_length(len_rounding_approx)
      )
      idx <- idx + 1L

      out[[idx]] <- data.frame(
        rep = b,
        p = p,
        method = "Rounding:corrected",
        covered = cov_rounding_corrected,
        length = clean_length(len_rounding_corrected)
      )
      idx <- idx + 1L

      # ------------------------------------------------------
      # Example 2: One-step update
      # ------------------------------------------------------
      cov_onestep_approx <- onegrad_contains_linear_abs(
        X, y, x_new, y_new, alpha, eta = eta_onegrad
      )
      len_onestep_approx <- onegrad_length_linear_abs(
        X, y, x_new, alpha, eta = eta_onegrad
      )

      cov_onestep_corrected <- corrected_onegrad_contains_linear_abs(
        X, y, x_new, y_new, alpha, eta = eta_corrected_onegrad
      )
      len_onestep_corrected <- corrected_onegrad_length_linear_abs(
        X, y, x_new, alpha, eta = eta_corrected_onegrad
      )

      out[[idx]] <- data.frame(
        rep = b,
        p = p,
        method = "One-step update:approx",
        covered = cov_onestep_approx,
        length = clean_length(len_onestep_approx)
      )
      idx <- idx + 1L

      out[[idx]] <- data.frame(
        rep = b,
        p = p,
        method = "One-step update:corrected",
        covered = cov_onestep_corrected,
        length = clean_length(len_onestep_corrected)
      )
      idx <- idx + 1L

      # ------------------------------------------------------
      # Example 3: Bayesian AOI-PPD
      # ------------------------------------------------------
      if (include_bayesian) {
        bayes_df <- bayesian_is_results_for_dataset(
          X = X,
          y = y,
          x_new = x_new,
          y_true = y_new,
          alpha = alpha,
          T = T_bayes,
          mu_prior = mu_prior,
          sig_prior = sig_prior,
          sig_lik_model = sig_lik_model,
          min_y_aoi = bayes_min_y_aoi,
          max_y_aoi = bayes_max_y_aoi,
          min_y_corrected = bayes_min_y_corrected,
          max_y_corrected = bayes_max_y_corrected,
          length_grid_size = bayes_length_grid_size,
          shrink = bayes_shrink
        )

        for (jj in seq_len(nrow(bayes_df))) {
          out[[idx]] <- data.frame(
            rep = b,
            p = p,
            method = bayes_df$method[jj],
            covered = bayes_df$covered[jj],
            length = clean_length(bayes_df$length[jj])
          )
          idx <- idx + 1L
        }
      }
    }
  }

  cat("\n")
  do.call(rbind, out[seq_len(idx - 1L)])
}

# ============================================================
# Summaries and plotting
# ============================================================

summarize_cp_results <- function(res) {
  res %>%
    group_by(p, method) %>%
    summarise(
      coverage = mean(as.numeric(covered), na.rm = TRUE),
      coverage_se = sd(as.numeric(covered), na.rm = TRUE) /
        sqrt(sum(!is.na(covered))),
      avg_length = mean(length, na.rm = TRUE),
      length_se = sd(length, na.rm = TRUE) /
        sqrt(sum(!is.na(length))),
      .groups = "drop"
    )
}

make_cp_family_plot <- function(summ, alpha = 0.10) {
  method_order <- c(
    "Deletion:approx", "Deletion:corrected",
    "Rounding:approx", "Rounding:corrected",
    "One-step update:approx", "One-step update:corrected",
    "Bayesian:approx", "Bayesian:corrected"
  )
  
  method_info <- data.frame(
    method = method_order,
    family = c(
      "Example 0: Deletion", "Example 0: Deletion",
      "Example 1: Rounding", "Example 1: Rounding",
      "Example 2: One-step update", "Example 2: One-step update",
      "Example 3: Bayesian", "Example 3: Bayesian"
    ),
    version = rep(c("approx", "corrected"), 4),
    stringsAsFactors = FALSE
  )
  
  summ_plot <- summ %>%
    left_join(method_info, by = "method") %>%
    filter(!is.na(family))
  
  plot_df <- bind_rows(
    summ_plot %>%
      transmute(
        p,
        method,
        family,
        version,
        metric = "Coverage",
        value = coverage,
        se = coverage_se
      ),
    summ_plot %>%
      transmute(
        p,
        method,
        family,
        version,
        metric = "Length",
        value = avg_length,
        se = length_se
      )
  )
  
  plot_df$family <- factor(
    plot_df$family,
    levels = c("Example 0: Deletion", "Example 1: Rounding", "Example 2: One-step update", "Example 3: Bayesian")
  )
  
  plot_df$metric <- factor(
    plot_df$metric,
    levels = c("Coverage", "Length")
  )
  
  plot_df$method <- factor(plot_df$method, levels = method_order)
  
  plot_df$version <- factor(
    plot_df$version,
    levels = c("approx", "corrected"),
    labels = c("Approximate", "Tournament")
  )
  
  version_cols <- c(
    "Approximate" = "#1f78b4",
    "Tournament" = "#e31a1c"
  )
  
  ref_df <- expand.grid(
    metric = factor("Coverage", levels = levels(plot_df$metric)),
    family = levels(plot_df$family)
  )
  ref_df$yint <- 1 - alpha
  
  ggplot(
    plot_df,
    aes(x = p, y = value, color = version, group = method)
  ) +
    geom_hline(
      data = ref_df,
      aes(yintercept = yint),
      linetype = "dashed",
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    geom_errorbar(
      aes(
        ymin = pmax(0, value - se),
        ymax = ifelse(
          metric == "Coverage",
          pmin(1, value + se),
          value + se
        )
      ),
      width = 0.8,
      alpha = 0.6
    ) +
    scale_color_manual(
      name = "Method",
      values = version_cols,
      breaks = c("Approximate", "Tournament")
    ) +
    facet_grid(metric ~ family, scales = "free_y", switch = "y") +
    labs(
      x = "p (dimension)",
      y = NULL
    ) +
    theme_bw(base_size = 26) +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey85"),
      strip.text = element_text(face = "bold", size = 24),
      axis.title = element_text(size = 26),
      axis.text = element_text(size = 26),
      legend.title = element_text(size = 26),
      legend.text = element_text(size = 26)
    )
}

# ============================================================
# Default run
# ============================================================

alpha <- 0.10
n <- 100
p_values <- seq(20, 100, 5)
reps <- 100
M <- 10
sigma <- 1
seed <- 1

eta_onegrad <- 10
eta_corrected_onegrad <- 10

T_bayes <- 100
mu_prior <- 10.0
sig_prior <- 1.0
sig_lik_model <- 1.0

res <- simulate_dim_curve_all(
  n = n,
  p_values = p_values,
  reps = reps,
  alpha = alpha,
  M = M,
  sigma = sigma,
  seed = seed,
  eta_onegrad = eta_onegrad,
  eta_corrected_onegrad = eta_corrected_onegrad,
  T_bayes = T_bayes,
  mu_prior = mu_prior,
  sig_prior = sig_prior,
  sig_lik_model = sig_lik_model,
  include_bayesian = TRUE
)

summ <- summarize_cp_results(res)
p_all <- make_cp_family_plot(summ, alpha = alpha)

print(p_all)

ggsave(
  file.path(script_dir, "corrected_vs_uncorrected_eta_10.pdf"),
  plot = p_all,
  width = 24,
  height = 8
)
