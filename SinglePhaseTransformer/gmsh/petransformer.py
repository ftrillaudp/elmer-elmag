### Frederic Trillaud <ftrillaudp@gmail.com>
### May 2025

from mpi4py import MPI
import gmsh
import sys
import numpy as np

name = "petransformer"

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
gdim = 2

flag_visu = 1

### Dimensions
corethickness = 0.099
cornerradius = 0.*0.002
corewidth = 0.065
coreheight = 0.067
coreaperturewidth = 0.012
coreapertureheight = 0.047
airgap = 0.002
corecentrallegwidth = 0.019
corelegwidth = 0.5*(corewidth-(corecentrallegwidth+2*coreaperturewidth))
airRadius = 0.2
shellRadius = 1.2*airRadius

strandRadius = 0.001
nbturns_1 = 9
distancecenterstrands = 2*strandRadius+0.001
layer1height = nbturns_1 * (2*strandRadius) + (nbturns_1-1) * distancecenterstrands
separationcorelayer_1 = 0.002
xposstrand_1 = 0.5*corecentrallegwidth+separationcorelayer_1+strandRadius

lc_1 = coreheight / 50
lc_2 = shellRadius / 10

### Get the tags for dimtags list
def get_gdimtags(dimtags, gdim):
	tags = list()
	for i, v in enumerate(dimtags):
		if (v[0] == gdim):
			tags.append(v[1])
	return tags	


if mesh_comm.rank == model_rank:
	corebulk = cascade.addRectangle(-0.5*corewidth, -0.5*coreheight, 0., corewidth, coreheight, 1, cornerradius)
	coretool_1 = cascade.addRectangle(-0.5*corecentrallegwidth-coreaperturewidth, -0.5*coreapertureheight, 0., coreaperturewidth, coreapertureheight, 2, cornerradius)
	coretool_2 = cascade.addRectangle(0.5*corecentrallegwidth, -0.5*coreapertureheight, 0., coreaperturewidth, coreapertureheight, 3, cornerradius)
	coretool_3 = cascade.addRectangle(-0.5*corewidth, -0.5*airgap, 0., corewidth, airgap, 4, 0.*cornerradius)
	core = cascade.cut([(gdim, 1)], [(gdim, i) for i in range(coretool_1,coretool_3+1)], removeObject = True, removeTool = True)
	core_tags = get_gdimtags(core[0], gdim)
	
	### Turns
	strand_1 = cascade.addDisk(xposstrand_1, 0., 0., strandRadius, strandRadius)
	print("strand_1", strand_1)
	
	print("core", core[0])
	
	airbulk = cascade.addDisk(0., 0., 0., airRadius, airRadius)
	air = cascade.cut([(gdim, airbulk)], core[0], removeObject = True, removeTool = False)
	air_tags = get_gdimtags(air[0], gdim)
	
	shellbulk = cascade.addDisk(0., 0., 0., shellRadius, shellRadius)
	shell = cascade.cut([(gdim, shellbulk)], [(gdim, i) for i in air_tags], removeObject = True, removeTool = False)
	shell_tags = get_gdimtags(shell[0], gdim)
	
	cascade.removeAllDuplicates
	
	cascade.synchronize()
	
	meshing.setSize(mod.getEntities(0), lc_1)
	shell_surfNodes = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=True)
	meshing.setSize(shell_surfNodes, lc_2)
	
	shell_surf = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=False)
	shell_dt = get_gdimtags(shell_surf, gdim-1)
	
	# ~ mod.addPhysicalGroup(gdim, core_tags, 1, name="core")
	# ~ mod.addPhysicalGroup(gdim, [7], 2, name="coilr")
	# ~ mod.addPhysicalGroup(gdim, [8], 3, name="coill")
	# ~ mod.addPhysicalGroup(gdim, air_tags, 6, name="air")
	# ~ mod.addPhysicalGroup(gdim, shell_tags[:-2], 7, name="shell")
	# ~ mod.addPhysicalGroup(gdim-1, [shell_dt[1]], 8, name="boundary")
	
	# ~ mod.setColor(core[0], 127, 127, 127, recursive=True)
	# ~ mod.setColor([(gdim, 7)], 255, 0, 0, recursive=False)
	# ~ mod.setColor([(gdim, 8)], 255, 51, 51, recursive=False)
	# ~ mod.setColor([(gdim, 9)], 204, 102, 0, recursive=False)
	# ~ mod.setColor([(gdim, 10)], 255, 153, 51, recursive=False)
	# ~ mod.setColor(air[0], 0, 128, 255, recursive=False)
	# ~ mod.setColor([(gdim, 16)], 153, 153, 255, recursive=False)
	# ~ mod.setColor([(gdim-1, shell_dt[1])], 255, 255, 255, recursive=False)
	
	meshing.removeDuplicateNodes
	meshing.removeDuplicateElements
	
	meshing.generate(gdim)
	
	gmsh.write(name+".msh")

gmsh.finalize()

# Build the input file for the solver	
# ~ myFile0 = open(name+".par", "w")
# ~ with myFile0 as file:
	# ~ file.write("coreID = 1;\n")
	# ~ file.write("coilID_HVr = 2;\n")
	# ~ file.write("coilID_HVl = 3;\n")
	# ~ file.write("coilID_LVr = 4;\n")
	# ~ file.write("coilID_LVl = 5;\n")
	# ~ file.write("airID = 6;\n")
	# ~ file.write("shellID = 7;\n")
	# ~ file.write("boundaryID = 8;\n")
	# ~ file.write("Ae_HV = "+str(hvwidth*hvheight)+";\n")
	# ~ file.write("Ae_LV = "+str(lvwidth*lvheight)+";\n")
	# ~ file.write("corethickness = "+str(corethickness)+";\n")
	# ~ file.write("airRadius = "+str(airRadius)+";\n"+"shellRadius = "+str(shellRadius)+";")
# ~ myFile0.close()

