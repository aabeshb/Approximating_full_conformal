## This file reads all the data from the results file which are obtained by running run.R
## and then creates the final summary table with the mean and standard error reported in 
## the table in the paper

## set dir_path to the directory containing all the result files produced by run.R

dir_path <- "~/Generalization of cross-conformal/Diabetes rounding/results"

files <- list.files(
  dir_path,
  pattern = "^res_[0-9]+\\.rds$",
  full.names = TRUE
)

length(files)   # should be 100
sort(basename(files))[1:5]

read_one <- function(f) {
  obj <- readRDS(f)
  
  # handle either a saved summary data frame directly,
  # or a saved list with a $summary element
  if (is.list(obj) && !is.data.frame(obj) && "summary" %in% names(obj)) {
    obj <- obj$summary
  }
  
  obj$source_file <- basename(f)
  obj
}

all_summaries <- do.call(rbind, lapply(files, read_one))

# optional: inspect combined data
print(head(all_summaries))
print(table(all_summaries$method))

# final aggregation across the 100 iterations
final_summary <- aggregate(
  cbind(mean_coverage, mean_length) ~ method,
  data = all_summaries,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE)/sqrt(length(x)))
)

# make it into a clean data frame
final_summary <- data.frame(
  method = final_summary$method,
  avg_coverage = final_summary$mean_coverage[, "mean"],
  se_coverage  = final_summary$mean_coverage[, "se"],
  avg_length   = final_summary$mean_length[, "mean"],
  se_length    = final_summary$mean_length[, "se"],
  row.names = NULL
)

print(final_summary)