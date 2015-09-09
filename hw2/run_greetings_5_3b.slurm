#!/bin/bash

#SBATCH --job-name=greetings
#SBATCH --output=bgreetings_5_3_slurm.out
#SBATCH --error=bgreetings_5_3_slurm.err
#SBATCH --partition=batch
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=3

srun ./greetings
