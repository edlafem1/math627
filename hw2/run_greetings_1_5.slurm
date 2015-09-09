#!/bin/bash

#SBATCH --job-name=greetings
#SBATCH --output=greetings_1_5_slurm.out
#SBATCH --error=greetings_1_5_slurm.err
#SBATCH --partition=develop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5

srun ./greetings
