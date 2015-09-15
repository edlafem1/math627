#!/bin/bash
#SBATCH --job-name=power
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=develop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --constraint=hpcf2013

srun ./power    4 1.0e-12 50

