# Analysis Code for the Bayesian Level-Set (BALLET) Clustering Article

This repository contains code necessary to reproduce results from the Bayesian Level-Set Clustering article by David Buch, Miheer Dewaskar, and David Dunson. 

The repository is organized into four top-level directories: `data`, `R`, `scripts`, and `src`. 

1. When cloned, the `data` folder will be nearly empty, containing only a copy of the t-SNE embedded RNA sequencing data found on [Renesh Bedre's website](https://www.reneshbedre.com/blog/dbscan-python.html) and the raw cosmological sky survey data, subsampled from the Edinburgh-Durham Southern Galaxy Catalogue, provided to us via personal correspondence from [Woncheol Jang](https://wcjang.github.io/). The various other datasets which appear in our article are simulated by the scripts in this repository, especially `scripts/data_processing/toy_data_processing.R`.

2. The top-level directory `R` contains `R` function definitions and other helper code that is used in the course of our BALLET data analyses.

3. Similar to the `R` directory, the `src` directory contains code for functions which are used to perform the data analyses. However, the functions defined in `src` are written in C++, and bound to `R` functions using the `Rcpp` package.

Together, the `R` and `src` directories contain a complete collection of functions needed to carry out a BALLET analysis of a new dataset, should a user be interested. We may use these two folders as the foundation for an open-source BALLET software package in  `R`.

4. The `scripts` directory contains the main `R` scripts which use the functions defined in `R` and `src` to carry out BALLET analyses and produce the output figures and tables seen in the main article and supplement. For ease of reading, the scripts in this top-level directory are grouped into four sub-folders: `data_processing`, `illustrations`, `toy_challenge_analysis`, and `sky_survey_analysis`.

