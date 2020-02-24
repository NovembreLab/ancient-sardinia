#!/bin/bash
#SBATCH --partition=jnovembre
#SBATCH --ntasks=1
#SBATCH --job-name="extractEigenStrat"
#SBATCH --time=2:00:00
#SBATCH --mem=56G
#SBATCH --mail-user=hringbauer@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=1

module load python/3.7.0
python cluster_run.py
