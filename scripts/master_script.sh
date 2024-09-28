#!/bin/bash

# --- Process relevant data files before running our analyses ---

Rscript ./scripts/data_processing/edsgc_processing.R
Rscript ./scripts/data_processing/toy_data_processing.R

# --- Generate the illustration in Figure 1 and a supplemntary toy figure.  ---

Rscript ./scripts/illustrations/gaussian_splitting_example.R
Rscript ./scripts/illustrations/univariate_different_levels.R

# --- Run the illustrative challenge datasets in Section 5 ---

## Plot the three challenge datasets

Rscript ./scripts/toy_challenge_analysis/toy_data.R

## Run our analysis with three density models: 
### 1. Dirichlet Process Mixture Model (DPMM),
### 2. Adaptive Polya Trees (APT), and 
### 3. NNDM (Nearest Neighbor Dirichlet Mixtures).

Rscript ./scripts/toy_challenge_analysis/toy_dpmm.R
Rscript ./scripts/toy_challenge_analysis/toy_apt.R
Rscript ./scripts/toy_challenge_analysis/toy_nndm.R

## Plot the results from clustering the challenge datasets.

Rscript ./scripts/toy_challenge_analysis/toy_plotting.R
 
# --- Analysis of the astronomical sky survey data in Section 6 ---

## This script estimates the parametes of DBSCAN to use for our analysis.

Rscript ./scripts/sky_survey_analysis/pre_simulation.R

## --- Analyze simulated data ---
### First, we run multiple instances of `single_simulation.R` in parallel using the 
####      [SLURM Workload Manager](https://slurm.schedmd.com/documentation.html).
####      see the file `run_simulations.sh` for an example.
### Then, we run the file `process_simulations.R` to gather results.

## Here is the equivalent place-holding code:

for i in {1..50}; do
  env SLURM_ARRAY_TASK_ID=$i Rscript ./scripts/sky_survey_analysis/single_simulation.R
done

Rscript ./scripts/sky_survey_analysis/process_simulations.R

## --- Analyze the Edinburgh-Durham Southern Galaxy Catalogue (EDSGC) ----

Rscript ./scripts/sky_survey_analysis/edsgc.R

