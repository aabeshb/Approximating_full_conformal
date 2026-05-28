# Approximating full conformal

Code for replicating experiments in the paper *"Approximating full conformal prediction: distribution
free guarantees via the tournament correction."*

## Experiments workflow instructions

This document provides instructions for running simulations and the real data experiments in the paper. For the simulations you can 
simply run the .R files and for the real data experiments you can use SLURM to submit jobs in batches. For users without SLURM, 
the corresponding `run.R` files can be executed in a simple loop with appropriate seed values.

---

### Simulation 1: Comparison of coverage and length

This experiment plots the coverage and the length of the prediction sets of both the approximate methods and the
tournament corrected methods in the paper.

1. Navigate to the directory `simulation_coverage_length`:
   ```bash
   cd simulation_coverage_length
   ```

2. Run the following R script to create the plot:
   ```bash
   Rscript plot.R
   ```

---

### Real data experiments 

This experiment produces the table in the paper that compares the performance of the approximate methods and the tournament methods
on the diabetes data for both the **Bayesian** and the **Rounding** method. 

1. Navigate to the appropriate directory (`Diabetes rounding`, or `Diabetes bayes`):
   ```bash
   cd Diabetes\ bayes
   # or: cd Diabetes\ rounding
   ```

2. Submit jobs using SLURM:
   ```bash
   sbatch --array=1-50 job.sbatch
   # or(rounding): --array=1-100 job.sbatch
   ```

3. Compile results:
   ```bash
   Rscript combine_results.R
   ```
---

### Simulation: stability comparison

This experiment compares the stability requirement for the approximate and the tournament methods.


1. Run the R script:
   ```bash
   Rscript stability_corrected_uncorrected.R
   ```
---



