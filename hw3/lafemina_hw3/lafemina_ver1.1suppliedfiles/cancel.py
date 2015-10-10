import os

file = open("all.txt", 'r')
for line in file:
    line = line.split()
    os.system("/cm/shared/apps/slurm/14.11.6/bin/scancel " + line[0])
file.close()


