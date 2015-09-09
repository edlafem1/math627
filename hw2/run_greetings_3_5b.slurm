#!/bin/bash

#SBATCH --job-name=greetings
#SBATCH --output=bgreetings_3_5_slurm.out
#SBATCH --error=bgreetings_3_5_slurm.err
#SBATCH --partition=batch
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=5

srun ./greetings
