#!/bin/bash

#SBATCH --job-name=qpAdm
#SBATCH --output=log/qpAdm-lm.%J.out
#SBATCH --error=log/qpAdm-lm.%J.err
#SBATCH --time=24:00:00
#SBATCH --partition=jnovembre
#SBATCH --nodes=1
#SBATCH --exclude=midway2-[0401-0402]
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2g

module load openblas
module load gsl

bin/qpAdm -p $1 > $2

