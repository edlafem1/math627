import sys

filename = "run_" + str(sys.argv[1]) + "_" + str(sys.argv[2]) + "_" + str(sys.argv[3]) + "_n=" + str(sys.argv[4]) + ".slurm"

slurm_file = open(filename, "w")

slurm_file.write("#!/bin/bash\n")
slurm_file.write("#SBATCH --job-name=" + str(sys.argv[1]) + str(sys.argv[4]) + "\n")
slurm_file.write("#SBATCH --output=" + str(sys.argv[1]) + "_" + str(sys.argv[2]) + "_" + str(sys.argv[3]) + "_n=" + str(sys.argv[4]) + ".out\n")
slurm_file.write("#SBATCH --error=" + str(sys.argv[1]) + "_" + str(sys.argv[2]) + "_" + str(sys.argv[3]) + "_n=" + str(sys.argv[4]) + ".err\n")
slurm_file.write("#SBATCH --partition=batch\n")
slurm_file.write("#SBATCH --nodes=" + str(sys.argv[2]) + "\n")
slurm_file.write("#SBATCH --ntasks-per-node=" + str(sys.argv[3]) + "\n")
slurm_file.write("#SBATCH --exclusive\n")
slurm_file.write("#SBATCH --constraint=hpcf2013\n")
slurm_file.write("#SBATCH --qos=short\n")
slurm_file.write("\n")
slurm_file.write("srun ./" + str(sys.argv[1]))
if len(sys.argv) > 4:
        for i in range(4,len(sys.argv)):
                slurm_file.write(" " + str(sys.argv[i]))

slurm_file.close()

import subprocess

sbatch = "/cm/shared/apps/slurm/14.11.6/bin/sbatch"
subprocess.call([sbatch, filename])
