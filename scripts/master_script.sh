#!/bin/bash


Rscript ./scripts/data_processing/edsgc_processing_wj.R
Rscript ./scripts/data_processing/toy_data_processing.R

Rscript ./scripts/illustrations/gaussian_splitting_example.R
Rscript ./scripts/illustrations/univariate_different_levels.R

Rscript ./scripts/toy_challenge_analysis/toy_data.R
Rscript ./scripts/toy_challenge_analysis/toy_dpmm.R
Rscript ./scripts/toy_challenge_analysis/toy_apt.R
Rscript ./scripts/toy_challenge_analysis/toy_nndm.R

Rscript ./scripts/toy_challenge_analysis/bound_analysis.R


