# F. Trillaud <ftrillaudp@gmail.com>
# 07/09/2024

### LIBRARIES ###
import time
import tracemalloc

tracemalloc.start()
st = time.time()

import os
import sys
import numpy as np
import random as rd
from termcolor import colored
import pathlib
import psutil
# ~ import re

st_1 = time.time()
print("Load libraries: {0:.6f} s".format(st_1-st))

### FUNCTION (UDF) ###
def varname(var, dir=locals()):
    return [key for key, val in dir.items() if id(val) == id(var)][0]

# Get the current directory path
myPath_case = str(pathlib.Path(__file__).parent.resolve())
myPath_resu = "./resu/"

os.system("rm -R "+myPath_resu+"; rm ./*.pre ./*.res ./*.par")

name = "transformer"
getdp = "/opt/onelab-Linux64/getdp "
cmd_msh = 'mpirun -n 1 python '+name+'.py -2 -bin -v 2 -nopopup'
cmd_solver = getdp+name+' -msh '+name+'.msh -bin -solve resolution -cpu -v 5'
cmd_post = getdp+name+' -msh '+name+'.msh -pos postOperation -cpu -v 5'

st_2 = time.time()
print("Initialization: {0:.6f} s".format(st_2-st_1))


### RUN ###
os.system(cmd_msh)
st_2 = time.time()
print("Geometry and meshing: {0:.6f} s".format(st_2-st_1))

os.system(cmd_solver)
os.system(cmd_post)
st_3 = time.time()
print("Runs: {0:.6f} s".format(st_3-st_2))


tracemalloc.stop()
process = psutil.Process(os.getpid())
et = time.time()
print("CPU time (s) = {0:.3f}".format(et-st))
print("Total memory (Mb) = {0:.3f}".format(1e-6*process.memory_info().rss))  # in bytes

