### Frederic Trillaud <ftrillaudp@gmail.com>
### May 2025

from mpi4py import MPI
import gmsh
import sys
import numpy as np

name = "coil"

model_rank = 0
mesh_comm = MPI.COMM_WORLD

gmsh.initialize(sys.argv)

mod = gmsh.model
cascade = mod.occ
meshing = mod.mesh
geom = mod.geo

mod.add(name)
mod.add("boolean")

gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.Algorithm3D', 10)

### Model dimension
gdim = 3

########################################
### Coil geometry and current densities
Ri = 0.05 # Inner coil radius of cylinder
wd = 0.02 # Coil width
Ro = Ri+wd # Outer radius of magnet
h = 0.02 # coil height

x = [0, 0]; dx = [0, 0]
y = [0, 0]; dy = [0, 0]
z = [0,0]; dz = [h, h]
r = [Ro, Ri]

### Mesh size for coil
cl = h/3;


### Get the tags for dimtags list
def get_gdimtags(dimtags, gdim):
	tags = list()
	for i, v in enumerate(dimtags):
		if (v[0] == gdim):
			tags.append(v[1])
	return tags		

### Main loop
if (mesh_comm.rank == model_rank):
    coil = list()
    for i,v in enumerate(x):
        coil.append(cascade.addCylinder(v, y[i], z[i], dx[i], dy[i], dz[i], r[i], i))
    cascade.cut([(3, 0)], [(3, 1)], 2, True, True)
    cascade.mesh.setSize(cascade.getEntities(0), cl)

    ### Air
    Ra = 2*Ro; lc_air = Ra / 8
    air = cascade.addSphere(0, 0, 0, Ra, 3)
    cascade.cut([(3, 3)], [(3, 2)], 4, True, False)

    ### Remove duplicate nodes
    cascade.removeAllDuplicates
    
    ### Sunchronize Opencascade engine with gmsh
    cascade.synchronize()
    
    ### Remove duplicate mesh nodes and elements
    meshing.removeDuplicateNodes
    meshing.removeDuplicateElements
    
    ### Generate mesh
    meshing.generate(gdim)
    
    ### Coil
    mod.addPhysicalGroup(gdim, [2], 1)
    mod.setColor([(gdim, 2)], 255, 165, 0)
    mod.setPhysicalName(gdim, 1, "coil")
    
    ### Air
    mod.addPhysicalGroup(gdim, [4], 2)
    mod.setColor([(gdim, 4)], 0, 191, 255)
    mod.setPhysicalName(gdim, 2, "air")

        
    ### Get the boundaries
    bnd_dt = gmsh.model.getBoundary(dimTags = [(gdim, 4)], combined = True, oriented = False, recursive = False)
    bnd_t = get_gdimtags(bnd_dt, gdim-1)
    print(bnd_dt[-1])
    ### Boundary id from the center of mass
    mod.addPhysicalGroup(gdim-1, [bnd_t[-1]], 3)
    mod.setColor([bnd_dt[-1]], 255, 0, 0)
    mod.setPhysicalName(gdim-1, 3, "airboundary")
    
    ### Save the mesh in gmsh format,. It will be read by the solver
    gmsh.write("mesh.msh")
    #gmsh.write(name+".step")

# ~ # Launch the GUI to see the results:
# ~ if '-nopopup' not in sys.argv:
    # ~ gmsh.fltk.run()

# Build the input file for the solver	
myFile0 = open(name+".dat", "w")
with myFile0 as file:
    file.write("$ Ri = "+str(Ri)+"\n")
    file.write("$ Ro = "+str(Ro)+"\n")
    file.write("$ wd = "+str(wd)+"\n")
    file.write("$ h = "+str(h))
myFile0.close()

### Close the gmsh engine	
gmsh.finalize()
