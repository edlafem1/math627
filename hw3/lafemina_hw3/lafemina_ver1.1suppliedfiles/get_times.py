
import os
from os import listdir
from os.path import isfile, join
mypath = "."
onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
'''
ns = [1024, 2048, 4096, 8192, 16384, 32768, 65536]
#ns = [65536]
total_p = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
ntasks_per_node = [1, 2, 4, 8, 8, 8, 8, 8, 8, 8, 16]
'''
n = dict()
err_file = open("error_list.txt", "w")

for filename in onlyfiles:
    if (os.path.splitext(filename)[1] == ".out"):
        str_filename = filename
        filename = filename.split("_")
        n_index = len(filename)-1
        filename[len(filename)-1] = filename[len(filename)-1][2:filename[len(filename)-1].find(".")]
        #print filename
        if filename[0] == "power" and len(filename) == 4 and int(filename[n_index]) >= 1024:
            if filename[n_index] not in n:
                n[filename[n_index]] = dict()
            p = int(filename[1]) * int(filename[2])
            res_file = open(str_filename, "r")
            lines = res_file.readlines()
            res_file.close()
            if (len(lines) == 0 or lines[-1][0] != "T"):
                err_file.write("Error: " + str_filename + "\n")
                continue
            time = lines[-1]
            time = time.split(" ")
            n[filename[n_index]][str(p)] = time[-1]
for val in n:
    print val
    print n[val]
    print "\n"
