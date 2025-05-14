### Frederic Trillaud <ftrillaudp@gmail.com>
### May 2025

from mpi4py import MPI
import gmsh
import sys
import numpy as np

gmsh.initialize()

name = "bulk"

model_rank = 0
mesh_comm = MPI.COMM_WORLD

gmsh.initialize(sys.argv)

mod = gmsh.model
cascade = mod.occ
meshing = mod.mesh
geom = mod.geo

mod.add(name)
mod.add("boolean")

# ~ gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
gmsh.option.setNumber('Mesh.Algorithm3D', 10)

### Model dimension
gdim = 3

flag_visu = 1

R = 15e-3
h = 15e-3
Ra = 0.1

lca = Ra/10
lcb = h/10


if mesh_comm.rank == model_rank:
	bulk = cascade.addCylinder(0., 0., -0.5*h, 0., 0., h, R, 1) # bulk
	cascade.addSphere(0., 0., 0., Ra, 2) # air
	
	air = cascade.cut([(gdim, 2)], [(gdim, 1)], 3, removeObject = True, removeTool = False)
	
	cascade.removeAllDuplicates
	
	cascade.synchronize()
	
	
	meshing.setSize(mod.getEntities(0), lca)
	bulk_surfNode_dimtag = mod.getBoundary([(gdim, 1)], combined=True, oriented=False, recursive=True)
	meshing.setSize(bulk_surfNode_dimtag, lcb)



	
	up, down = mod.getAdjacencies(gdim, 3)
	
	mod.addPhysicalGroup(gdim, [1], 1, name="bulk")
	mod.setColor([(gdim, 1)], 127, 127, 127, recursive=True)
	mod.addPhysicalGroup(gdim, [3], 2, name="air")
	mod.setColor([(gdim, 3)], 0, 0, 255, recursive=True)
	mod.addPhysicalGroup(gdim-1, [down[1]], 3, name="bulkBoundary")
	mod.setColor([(gdim-1, down[1])], 127, 127, 127, recursive=True)
	mod.addPhysicalGroup(gdim-1, [down[0]], 4, name="airBoundary")
	mod.setColor([(gdim-1, down[0])], 0, 0, 255, recursive=True)
	
	# ~ meshing.removeDuplicateNodes
	# ~ meshing.removeDuplicateElements
	
	meshing.generate(gdim)
	 
	gmsh.write("mesh.msh")

gmsh.finalize()

