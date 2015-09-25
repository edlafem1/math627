import sys

ns = [1024, 2048, 4096, 8192, 16384, 32768, 65536]
ns = [65536]
total_p = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
ntasks_per_node = [1, 2, 4, 8, 8, 8, 8, 8, 8, 8, 16]

for n in ns:
    for index_p in range(len(total_p)):
        p = total_p[index_p]
        ntasks = ntasks_per_node[index_p]
        nodes = p / ntasks

        filename = "run_" + str(sys.argv[1]) + "_" + str(nodes) + "_" + str(ntasks) + "_n=" + str(n) + ".slurm"

        slurm_file = open(filename, "w")

        slurm_file.write("#!/bin/bash\n")
        slurm_file.write("#SBATCH --job-name=" + str(sys.argv[1]) + str(n) + "\n")
        slurm_file.write("#SBATCH --output=" + str(sys.argv[1]) + "_" + str(nodes) + "_" + str(ntasks) + "_n=" + str(n) + ".out\n")
        slurm_file.write("#SBATCH --error=" + str(sys.argv[1]) + "_" + str(nodes) + "_" + str(ntasks) + "_n=" + str(n) + ".err\n")
        slurm_file.write("#SBATCH --partition=batch\n")
        slurm_file.write("#SBATCH --nodes=" + str(nodes) + "\n")
        slurm_file.write("#SBATCH --ntasks-per-node=" + str(ntasks) + "\n")
        slurm_file.write("#SBATCH --exclusive\n")
        slurm_file.write("#SBATCH --constraint=hpcf2013\n")
        slurm_file.write("#SBATCH --qos=short\n")
        slurm_file.write("#SBATCH --mem=63000\n")
        slurm_file.write("\n")
        slurm_file.write("srun ./" + str(sys.argv[1]) + " " + str(n) + " 1e-12 50")
        slurm_file.write("\n")

        slurm_file.close()

        import subprocess

        sbatch = "/cm/shared/apps/slurm/14.11.6/bin/sbatch"
        subprocess.call([sbatch, filename])
