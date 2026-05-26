# ============================================================
# Run one Diabetes train/test split for a single Slurm-array seed
# using:
#   1) grid-based Full CP with Huber-RLM
#   2) tournament-corrected round-one-at-a-time with Huber-RLM
#
# Intended usage:
#   sbatch --array=1-100 job.sbatch
# where each array task sets SLURM_ARRAY_TASK_ID to one seed.
#

# ============================================================

seed = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (!is.finite(seed) || is.na(seed)) {
  stop("SLURM_ARRAY_TASK_ID is not set. Run with sbatch --array=1-50 job.sbatch, or export SLURM_ARRAY_TASK_ID manually.", call. = FALSE)
}
seed <- as.integer(seed)
#seed <- 1
# Locate the actual directory containing the run.R file.
#

get_script_path <- function(target_file = "run.R", required_neighbor = "diabetes_rounding.R") {
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
  
  # 1. Explicit overrides. These must be first because interactive R/RStudio
  # sessions often have no reliable script path.
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

source(file.path(script_dir, "diabetes_rounding.R"))

# ---------- settings ----------
# One split per Slurm-array task; use m = 100 test points.
n_iter <- 1L
test_size <- 100L
alpha <- 0.1
thres <- 1
lamb <- 2
grid_size <- 100L
maxit <- 1000L
reltol <- 1e-9
verbose <- TRUE
print_progress <- TRUE
progress_every <- 1L

# ---------- run one split ----------
message("Starting seed ", seed, " with m = ", test_size, " test points")
message("Progress printing is ", if (isTRUE(print_progress)) "ON" else "OFF", "; progress_every = ", progress_every)
results <- run_diabetes_experiment(
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

# ---------- add seed to outputs ----------
results$summary$seed <- seed

# reorder columns so seed appears first when possible
results$summary <- results$summary[, c("seed", setdiff(names(results$summary), "seed")), drop = FALSE]

# ---------- save outputs ----------
out_dir <- file.path(script_dir, "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

seed_tag <- sprintf("%03d", seed)

out_file <- file.path(out_dir, sprintf("res_%d.rds", seed))
saveRDS(results$summary, out_file)


#print(results$summary)
