#!/bin/bash

#SBATCH --job-name=greetings
#SBATCH --output=greetings_3_5_slurm.out
#SBATCH --error=greetings_3_5_slurm.err
#SBATCH --partition=develop
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=5

srun ./greetings
