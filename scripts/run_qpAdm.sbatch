#!/bin/bash
#SBATCH --job-name=qpAdm
#SBATCH --output=log/qpAdm.%J.out
#SBATCH --error=log/qpAdm.%J.err
#SBATCH --time=1:00:00
#SBATCH --partition=jnovembre
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6g

module load openblas
module load gsl

bin/qpAdm -p $1 > $2
