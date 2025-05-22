### Frederic Trillaud <ftrillaudp@gmail.com>
### May 2025

from mpi4py import MPI
import gmsh
import sys
import numpy as np

name = "transformer"

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
corethickness = 0.099;
cornerradius = 0.0025
corewidth = 0.23
coreheight = 0.204
coreaperturewidth = 0.051
coreapertureheight = 0.14
corecentrallegwidth = 0.064
corelegwidth = 0.5*(corewidth-(corecentrallegwidth+2*coreaperturewidth))
hvwidth = 0.0174
hvheight = 0.1216
separationcorehv = 0.0022
lvwidth = 0.0121
lvheight = 0.1085
separationcorelv = 0.0022
airRadius = 0.5
shellRadius = 1.2*airRadius

xposhole = 0.5*(corecentrallegwidth+corelegwidth)+coreaperturewidth
yposhole = 0.5*(coreapertureheight+corelegwidth)
holeRadius = 0.006

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
	coretool_3 = cascade.addDisk(xposhole, yposhole, 0., holeRadius, holeRadius, 4)
	coretool_4 = cascade.addDisk(-xposhole, yposhole, 0., holeRadius, holeRadius, 5)
	core = cascade.cut([(gdim, 1)], [(gdim, i) for i in range(2,6)], 6, removeObject = True, removeTool = False)
	
	hv = cascade.addRectangle(0.5*corecentrallegwidth+separationcorehv, -0.5*hvheight, 0., hvwidth, hvheight, 7, cornerradius/2)
	hv = cascade.addRectangle(-0.5*corecentrallegwidth-separationcorehv-hvwidth, -0.5*hvheight, 0., hvwidth, hvheight, 8, cornerradius/2)
	
	lv = cascade.addRectangle(0.5*corecentrallegwidth+separationcorehv+hvwidth+separationcorelv, -0.5*lvheight, 0., lvwidth, lvheight, 9, cornerradius/2)
	lv = cascade.addRectangle(-0.5*corecentrallegwidth-separationcorehv-hvwidth-separationcorelv-lvwidth, -0.5*lvheight, 0., lvwidth, lvheight, 10, cornerradius/2)
	
	airbulk = cascade.addDisk(0., 0., 0., airRadius, airRadius)
	air = cascade.cut([(gdim, airbulk)], [(gdim, i) for i in range(6,11)], removeObject = True, removeTool = False)
	
	shellbulk = cascade.addDisk(0., 0., 0., shellRadius, shellRadius)
	shell = cascade.cut([(gdim, shellbulk)], [(gdim, i) for i in range(6,15)], removeObject = True, removeTool = False)
	
	cascade.removeAllDuplicates
	
	cascade.synchronize()
	
	meshing.setSize(mod.getEntities(0), lc_1)
	shell_surfNodes = mod.getBoundary([(gdim, 16)], combined=True, oriented=False, recursive=True)
	meshing.setSize(shell_surfNodes, lc_2)
	
	shell_surf = mod.getBoundary([(gdim, 16)], combined=True, oriented=False, recursive=False)
	shell_dt = get_gdimtags(shell_surf, gdim-1)
	
	mod.addPhysicalGroup(gdim, [6], 1, name="core")
	mod.addPhysicalGroup(gdim, [7], 2, name="HVr")
	mod.addPhysicalGroup(gdim, [8], 3, name="HVl")
	mod.addPhysicalGroup(gdim, [9], 4, name="HLr")
	mod.addPhysicalGroup(gdim, [10], 5, name="HLl")
	mod.addPhysicalGroup(gdim, [i for i in range(11,16)], 6, name="air")
	mod.addPhysicalGroup(gdim, [16], 7, name="shell")
	mod.addPhysicalGroup(gdim-1, [shell_dt[1]], 8, name="boundary")
	
	mod.setColor([(gdim, 6)], 127, 127, 127, recursive=False)
	mod.setColor([(gdim, 7)], 255, 0, 0, recursive=False)
	mod.setColor([(gdim, 8)], 255, 51, 51, recursive=False)
	mod.setColor([(gdim, 9)], 204, 102, 0, recursive=False)
	mod.setColor([(gdim, 10)], 255, 153, 51, recursive=False)
	mod.setColor(air[0], 0, 128, 255, recursive=False)
	mod.setColor([(gdim, 16)], 153, 153, 255, recursive=False)
	mod.setColor([(gdim-1, shell_dt[1])], 255, 255, 255, recursive=False)
	
	meshing.removeDuplicateNodes
	meshing.removeDuplicateElements
	
	meshing.generate(gdim)
	
	gmsh.write("transformer.msh")

gmsh.finalize()

# Build the input file for the solver	
myFile0 = open(name+".par", "w")
with myFile0 as file:
	file.write("coreID = 1;\n")
	file.write("coilID_HVr = 2;\n")
	file.write("coilID_HVl = 3;\n")
	file.write("coilID_LVr = 4;\n")
	file.write("coilID_LVl = 5;\n")
	file.write("airID = 6;\n")
	file.write("shellID = 7;\n")
	file.write("boundaryID = 8;\n")
	file.write("Ae_HV = "+str(hvwidth*hvheight)+";\n")
	file.write("Ae_LV = "+str(lvwidth*lvheight)+";\n")
	file.write("corethickness = "+str(corethickness)+";\n")
	file.write("airRadius = "+str(airRadius)+";\n"+"shellRadius = "+str(shellRadius)+";")
myFile0.close()

