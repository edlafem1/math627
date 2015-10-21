ppns = [1, 2, 2, 2, 2, 2]
cpus = [1, 1, 2, 4, 8, 16]

for i in range(len(cpus)):
    for mkn in [1024, 2048, 4096, 8192, 16384, 32768]:
        filename = str(mkn) + "run_power_" + str(cpus[i]) + "_" + str(ppns[i]) + ".slurm"
        slurm = open(filename, 'w')
        slurm.write("#!/bin/bash\n")
        slurm.write("#SBATCH --job-name=p" + str(mkn) + "\n")
        slurm.write("#SBATCH --output=results/" + str(mkn) + "power_" + str(cpus[i]) + "_" + str(ppns[i]) + ".out\n")
        slurm.write("#SBATCH --error=results/" + str(mkn) + "power_" + str(cpus[i]) + "_" + str(ppns[i]) + ".err\n")
        slurm.write("#SBATCH --partition=batch\n")
        slurm.write("#SBATCH --nodes=" + str(cpus[i]) + "\n")
        slurm.write("#SBATCH --ntasks-per-node=" + str(ppns[i]) + "\n")
        slurm.write("#SBSTCH --exclusive\n")
        slurm.write("#SBATCH --constraint=hpcf2013\n")
        slurm.write("#SBATCH --qos=medium\n")
        slurm.write("#SBATCH --mem=63000\n")
        slurm.write("srun ./power " + str(mkn) + " " + str(mkn) + " " + str(mkn) + "\n")
        slurm.close()
    
        import subprocess

        sbatch = "/cm/shared/apps/slurm/14.11.6/bin/sbatch"
        subprocess.call([sbatch, filename])
