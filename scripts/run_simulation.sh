#!/bin/bash
#SBATCH -o /hpc/group/dunsonlab/md440/logs/BALLET_sim_%A_%a.out
#SBATCH -e /hpc/group/dunsonlab/md440/logs/BALLET_sim_%A_%a.err
#SBATCH --job-name=BALLET
#SBATCH --mail-user=miheer.dewaskar@duke.edu
#SBATCH --mail-type=END,FAIL  
#SBATCH --array=1-50
#SBATCH --mem=40G  
#SBATCH --partition dunsonlab
#SBATCH --account dunsonlab

echo "Started $SLURM_ARRAY_TASK_ID"
cd /hpc/group/dunsonlab/md440/Projects/ballet_article
module load R
Rscript -e "source('scripts/sky_survey_analysis/single_simulation.R', echo=TRUE)" 
echo "Finished $SLURM_ARRAY_TASK_ID"
# See https://hpc.nmsu.edu/discovery/slurm/job-arrays/ for 
# information on running parallel scripts in R.


