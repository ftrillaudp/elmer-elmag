# F. Trillaud <ftrillaudp@gmail.com>
# 05/04/2022


### LIBRARIES ###
import os
import sys
import pathlib
import time


#filename = input("Enter the name of the Gmsh mesh file without extension: ")

# Get the current directory path
myPath = str(pathlib.Path(__file__).parent.resolve())

# SOLVERS #
#pathname, extension = os.path.splitext(sys.argv[1])
#filenameList = pathname.split('.')
#filename = filenameList[-1]
filename = "mesh"

# print("Path name:", pathname)
# print("extension:", extension)
# print("File name:", filename)

# Commands
cmd_msh = "ElmerGrid 8 2 "+filename+".unv -out mesh"
cmd_solver = "ElmerSolver case.sif"

st = time.time()

print("Mesh generation")
os.system(cmd_msh)

print("Solver")
os.system(cmd_solver)

et = time.time()
print("CPU time (s) =", et-st)
