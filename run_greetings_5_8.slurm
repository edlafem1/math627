#!/bin/bash

#SBATCH --job-name=greetings
#SBATCH --output=greetings_5_8_any_slurm.out
#SBATCH --error=greetings_5_8_any_slurm.err
#SBATCH --partition=develop
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=8

srun ./greetings
