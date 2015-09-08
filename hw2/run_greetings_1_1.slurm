#!/bin/bash

#SBATCH --job-name=greetings
#SBATCH --output=greetings_1_1_slurm.out
#SBATCH --error=greetings_1_1_slurm.err
#SBATCH --partition=develop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

srun ./greetings
