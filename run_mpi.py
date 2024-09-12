import os
import sys

num_processes = 4
code_name = "shear_evp"
# code_name = "droplet_ndi_Oh"

if( code_name.endswith(".c") ):
    code_name = code_name[:-2]

if(os.path.isfile(code_name)):
    os.remove(code_name)
if(os.path.isfile("_%s.c" % (code_name))):
    os.remove("_%s.c" % (code_name))
    
# os.system("qcc -source -D_MPI=1 %s.c" % (code_name))
# os.system("mpicc -Wall -std=c99 -O2 -D_XOPEN_SOURCE=700 _%s.c -o %s -lm"  % (code_name, code_name))
# string_command = ("mpiexec -np %d ./%s " % (num_processes, code_name))+ " ".join(sys.argv[1:])
# os.system(string_command)

os.system("qcc -fopenmp -O2 -Wall %s.c -lm -o %s" % (code_name, code_name))
os.system("export OMP_NUM_THREADS=%d" % (num_processes))
os.system("./%s " % (code_name) + " ".join(sys.argv[1:]) )



# mpiexec -np 2 ./filaments 0 0.1 5 0 1 100 0.000001 0.001 0.0001 0.0175 10 
