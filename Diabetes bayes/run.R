seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.finite(seed) || is.na(seed)) {
  stop("SLURM_ARRAY_TASK_ID is not set. Run with sbatch --array=1-50 job.sbatch, or export SLURM_ARRAY_TASK_ID manually.", call. = FALSE)
}
seed <- as.integer(seed)
#seed <- 1
# Locate the actual directory containing this run.R file.

get_script_path <- function(target_file = "run.R", required_neighbor = "diabetes_source.R") {
  to_abs <- function(path) {
    path <- path[1L]
    if (is.na(path) || !nzchar(path)) return(NA_character_)
    path <- path.expand(path)
    is_abs <- grepl("^/", path) || grepl("^[A-Za-z]:[\\/]", path)
    if (!is_abs) path <- file.path(getwd(), path)
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }

  valid_script <- function(path) {
    path <- to_abs(path)
    if (is.na(path) || !nzchar(path)) return(FALSE)
    file.exists(path) &&
      basename(path) == target_file &&
      file.exists(file.path(dirname(path), required_neighbor))
  }

  candidates <- character(0)

  # 1. Explicit overrides. 
  env_path <- Sys.getenv("RUN_R_PATH", unset = "")
  if (nzchar(env_path)) candidates <- c(candidates, env_path)

  env_dir <- Sys.getenv("RUN_R_DIR", unset = "")
  if (nzchar(env_dir)) candidates <- c(candidates, file.path(env_dir, target_file))

  # A more generic name, useful outside SLURM.
  project_dir <- Sys.getenv("PROJECT_DIR", unset = "")
  if (nzchar(project_dir)) candidates <- c(candidates, file.path(project_dir, target_file))

  # 2. Standard case: Rscript /absolute/path/to/run.R.
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    candidates <- c(candidates, sub("^--file=", "", file_arg[length(file_arg)]))
  }

  # 3. RStudio active source file. This is useful if you are pressing Run
  # from an open run.R file rather than using Rscript.
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- try(rstudioapi::getSourceEditorContext(), silent = TRUE)
    if (!inherits(ctx, "try-error") && !is.null(ctx$path) && nzchar(ctx$path)) {
      candidates <- c(candidates, ctx$path)
    }
  }

  # 4. If run via source('/path/to/run.R'), R often stores the sourced file in
  # the call stack as `ofile`. This only works while the file is being sourced.
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) candidates <- c(candidates, ofile)
  }

  # 5. SLURM and working-directory fallbacks. These are intentionally late
  # because they can point to the submission directory or an unrelated project.
  slurm_submit <- Sys.getenv("SLURM_SUBMIT_DIR", unset = "")
  if (nzchar(slurm_submit)) candidates <- c(candidates, file.path(slurm_submit, target_file))

  candidates <- c(candidates, file.path(getwd(), target_file))
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  for (path in candidates) {
    if (valid_script(path)) return(to_abs(path))
  }

  msg <- paste(
    "Could not locate the directory containing run.R and diabetes_source.R.",
    "Tried candidates:",
    paste(sprintf("  - %s", candidates), collapse = "\n"),
    "",
    "Fix one of the following:",
    "  1. Run from the actual folder: cd '/path/to/Diabetes_cache'; Rscript run.R",
    "  2. Run with absolute path: Rscript '/path/to/Diabetes_cache/run.R'",
    "  3. Or set an override before running:",
    "       export RUN_R_DIR='/path/to/Diabetes_cache'",
    "       export RUN_R_PATH='/path/to/Diabetes_cache/run.R'",
    sep = "\n"
  )
  stop(msg, call. = FALSE)
}

get_script_dir <- function() dirname(get_script_path())

make_script_relative_path <- function(path, script_dir) {
  path <- path.expand(path)
  if (grepl("^/", path) || grepl("^[A-Za-z]:[\\\\/]", path)) return(path)
  file.path(script_dir, path)
}

script_path <- get_script_path()
script_dir <- dirname(script_path)

# Keep execution sequential even when chains = 4.
options(mc.cores = 1)

source(file.path(script_dir, "diabetes_source.R"))
rstan::rstan_options(auto_write = TRUE)

dat <- load_diabetes_lars()
sp <- make_train_test_split(dat$X, dat$y, train_frac = 0.7, seed = seed)

results_dir <- make_script_relative_path(Sys.getenv("RESULTS_DIR", unset = "results"), script_dir)
cache_root <- make_script_relative_path(Sys.getenv("CACHE_DIR", unset = "cache"), script_dir)
cache_dir <- file.path(cache_root, sprintf("seed_%03d", seed))

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

start = Sys.time()

res <- run_bayes_lasso_conformal_split(
  seed = seed,
  X_train = sp$X_train,
  y_train = sp$y_train,
  X_test = sp$X_test,
  y_test = sp$y_test,
  alpha = 0.20,
  iter = 3000,
  warmup = 1000,
  chains = 4,
  loo_iter = 3000,
  loo_warmup = 1000,
  loo_chains = 4,
  grid_size = 100,
  max_test = NULL,
  cache_dir = cache_dir,
  corrected = TRUE,
  refresh = 0,
  verbose = FALSE,
  print_progress = TRUE
)
end <- Sys.time()
#print(end - start)
out_file <- file.path(results_dir, sprintf("res_%d.rds", seed))
saveRDS(res, out_file)

# print(summarize_bayes_results(res))
