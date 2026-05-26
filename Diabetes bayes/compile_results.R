## This file reads all the data from the results file which are obtained by running run.R
## and then creates the final summary table with the mean and standard error reported in 
## the table in the paper

## set dir_path to the directory containing all the result files produced by run.R

dir_path <- "~/Generalization of cross-conformal/Diabetes bayes/results"

files <- list.files(
  path = dir_path,
  pattern = "^res_[0-9]+\\.rds$",
  full.names = TRUE
)

# Compute coverage and length average within each seed/job
per_seed <- do.call(rbind, lapply(files, function(f) {
  res <- readRDS(f)
  seed <- as.numeric(gsub("res_|\\.rds", "", basename(f)))
  
  out <- aggregate(
    cbind(
      coverage = as.numeric(covered),
      length = length
    ) ~ method,
    data = res,
    FUN = mean
  )
  
  out$seed <- seed
  out
}))

# Compute mean and standard errors across the seed-level averages
final_summary <- do.call(rbind, lapply(split(per_seed, per_seed$method), function(d) {
  data.frame(
    method = unique(d$method),
    mean_coverage = mean(d$coverage),
    se_coverage = sd(d$coverage) / sqrt(nrow(d)),
    mean_length = mean(d$length),
    se_length = sd(d$length) / sqrt(nrow(d)),
    n_seeds = nrow(d)
  )
}))

final_summary <- final_summary[order(final_summary$method), ]

print(per_seed)
print(final_summary)
