#!/bin/bash

# See https://www.bolzanor.eu/post/2021-02-05-parallel_cmd/2021-02-05-parallel_cmd/
# for information on running parallel scripts in R.

mkdir -p output/sky_survey_analysis/sim_study/log/

for id in {1..25}
do
  #echo nohup_$id.out
  nohup Rscript --verbose scripts/sky_survey_analysis/single_simulation.R --args $id  &> output/sky_survey_analysis/sim_study/log/output_$id.out &
done
