### Frederic Trillaud <ftrillaudp@gmail.com>
### September 2025

from mpi4py import MPI
import gmsh
import sys
import numpy as np

name = "wire"

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
gmsh.option.setColor("General.BackgroundGradient", 255, 255, 255)
# ~ gmsh.option.setColor("General.Background", 255, 255, 255)
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
# ~ gmsh.option.setNumber('Mesh.Algorithm3D', 10)
gmsh.option.setNumber('Mesh.OptimizeThreshold', 0.3)

### Model dimension
gdim = 2

flag_visu = 1

### Dimensions
Rstrand = 2.5e-3
Rturn = 0.02
Xcore = 0.5*Rturn
Ycore = 0.5*Rturn
Zcore = 0.01
Rair = 2*max(Rturn, np.sqrt(Xcore**2+Ycore**2))
Rshell = 1.2*Rair
Rcorner = 2e-3
Rfillet = 2e-3

# Number of turns
nbw = 1

### Mesh sizes, defined on nodes
lc_all = Ycore / 50
lc_shell = Rshell / 10
lc_strands = Rstrand / 10
lc_core = Xcore / 10

### Get the tags for dimtags list
def get_gdimtags(dimtags, gdim):
	tags = list()
	for i, v in enumerate(dimtags):
		if (v[0] == gdim):
			tags.append(v[1])
	return tags	


if mesh_comm.rank == model_rank:
	### Core
	core_tag = cascade.addRectangle(-0.5*Xcore, -0.5*Ycore, 0., Xcore, Ycore, 1, Rcorner)
	
	# ~ ### Turns
	strand_p_tag = cascade.addDisk(Rturn, 0., 0., Rstrand, Rstrand, core_tag+1)
	strand_n_tag = cascade.addDisk(-Rturn, 0., 0., Rstrand, Rstrand, strand_p_tag+1)
	strands_tags = [strand_p_tag, strand_n_tag]
	
	toollist = [(gdim, core_tag), (gdim, strand_p_tag), (gdim, strand_n_tag)]
	
	### Air
	airbulk = cascade.addDisk(0., 0., 0., Rair, Rair)
	air_dimtag = cascade.cut([(gdim, airbulk)], toollist, removeObject = True, removeTool = False)
	air_tag = get_gdimtags(air_dimtag[0], gdim)[0]
	
	toollist.append(air_dimtag[0][0])
	
	### Shell
	shellbulk = cascade.addDisk(0., 0., 0., Rshell, Rshell)
	shell_dimtag = cascade.cut([(gdim, shellbulk)], toollist, removeObject = True, removeTool = False)
	shell_tag = get_gdimtags(shell_dimtag[0], gdim)[0]
	
	### Remove duplicates
	cascade.removeAllDuplicates
	
	### Symchronize OCCT with Gmsh
	cascade.synchronize()

	### Boundaries (gdim-1)
	shell_bnd_dimtag = mod.getBoundary(shell_dimtag[0], combined=True, oriented=False, recursive=False)
	shell_bnd_tag = get_gdimtags(shell_bnd_dimtag, gdim-1)[1]
	
	strands_bnd_dimtag = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_tags)], combined=True, oriented=False, recursive=False)
	strands_bnd_tags = get_gdimtags(strands_bnd_dimtag, gdim-1)

	### Mesh sizes (gdim-2)
	meshing.setSize(mod.getEntities(0), lc_all)
	
	shell_bndnodes = mod.getBoundary(shell_dimtag[0], combined=True, oriented=False, recursive=True)
	meshing.setSize(shell_bndnodes, lc_shell)
	
	strands_bndnodes = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_tags)], combined=True, oriented=False, recursive=True)
	meshing.setSize(strands_bndnodes, lc_strands)
	
	core_bndnodes = mod.getBoundary([(gdim, core_tag)], combined=True, oriented=False, recursive=True)
	meshing.setSize(core_bndnodes, lc_core)


	### Remove duplicate nodes and elements
	meshing.removeDuplicateNodes
	meshing.removeDuplicateElements

	### Define group of geometries and ID
	k = 1
	mod.addPhysicalGroup(gdim, [core_tag], k, name="core")
	k = k+1
	mod.addPhysicalGroup(gdim, [air_tag], k, name="air")
	k = k+1
	mod.addPhysicalGroup(gdim, [shell_tag], k, name="shell")
	k = k+1
	mod.addPhysicalGroup(gdim-1, [shell_bnd_tag], k, name="boundary")
	k = k+1
	mod.addPhysicalGroup(gdim, [strands_tags[0]], k, name="wirep")
	k = k+1
	mod.addPhysicalGroup(gdim, [strands_tags[1]], k, name="wiren")
	
	### Colors
	mod.setColor([(gdim, core_tag)], 60, 60, 60, recursive=False)
	mod.setColor([(gdim, air_tag)], 50, 50, 255, recursive=False)
	mod.setColor(shell_dimtag[0], 0, 100, 200, recursive=False)
	mod.setColor(shell_bnd_dimtag, 150, 0, 75, recursive=False)
	mod.setColor([(gdim, v) for i, v in enumerate(strands_tags)], 255, 0, 0, recursive=False)
	
	meshing.generate(gdim)

	gmsh.write(name+".msh")

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

# ~ # Build the input file for the solver	
myFile0 = open(name+".par", "w")
with myFile0 as file:
	file.write("airRadius = "+str(Rair)+";\n"+"shellRadius = "+str(Rshell)+";\n")
	file.write("nbw = "+str(nbw)+";\n")
	file.write("nbID = "+str(k)+";\n")
	file.write("corethickness = 0.01;\n")
	file.write("coreID = 1;\n")
	file.write("airID = 2;\n")
	file.write("shellID = 3;\n")
	file.write("boundaryID = 4;\n")
	file.write("wirepID_1 = 5;\n")
	file.write("wirenID_1 = 6;\n")
myFile0.close()

# ~ ### For elmerfem if necessary
# ~ myFile1 = open(name+".dat", "w")
# ~ with myFile1 as file:
	# ~ file.write("nbw = "+str(nbw)+";\n")
	# ~ file.write("nbID = "+str(k)+";\n")
# ~ myFile2.close()

