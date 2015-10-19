

for mkn in [1024, 2048, 4096, 8192, 16384, 32768]:
    filename = "run_power_1_1_" + str(mkn) + ".slurm"
    slurm = open(filename, 'w')
    slurm.write("#!/bin/bash\n")
    slurm.write("#SBATCH --job-name=p" + str(mkn) + "\n")
    slurm.write("#SBATCH --output=results/" + str(mkn) + "power_1_1.out\n")
    slurm.write("#SBATCH --error=results/" + str(mkn) + "power_1_1.err\n")
    slurm.write("#SBATCH --partition=batch\n")
    slurm.write("#SBATCH --nodes=1\n")
    slurm.write("#SBATCH --ntasks-per-node=1\n")
    slurm.write("#SBSTCH --exclusive\n")
    slurm.write("#SBATCH --constraint=hpcf2013\n")
    slurm.write("#SBATCH --qos=medium\n")
    slurm.write("#SBATCH --mem=63000\n")
    slurm.write("srun ./power " + str(mkn) + " " + str(mkn) + " " + str(mkn) + "\n")
    slurm.close()
    
    import subprocess

    sbatch = "/cm/shared/apps/slurm/14.11.6/bin/sbatch"
    subprocess.call([sbatch, filename])
