### Frederic Trillaud <ftrillaudp@gmail.com>
### May 2025

from mpi4py import MPI
import gmsh
import sys
import numpy as np
from doublehelix import double_helix

name = "3dpetransformer"

model_rank = 0
mesh_comm = MPI.COMM_WORLD

gmsh.initialize(sys.argv)

opt = gmsh.option
mod = gmsh.model
occ = mod.occ
meshing = mod.mesh
geom = mod.geo

mod.add(name)
mod.add("boolean")

# ~ gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
opt.setColor("General.BackgroundGradient", 255, 255, 255)
# ~ gmsh.option.setColor("General.Background", 255, 255, 255)
opt.setNumber('Mesh.Optimize', 1)
opt.setNumber('Mesh.OptimizeNetgen', 1)
opt.setNumber('Mesh.Algorithm3D', 10)
opt.setNumber('Mesh.OptimizeThreshold', 0.4)
# ~ opt.setNumber('Mesh.AngleToleranceFacetOverlap', 0.1)


### Get the tags for dimtags list
def get_gdimtags(dimtags, gdim):
	tags = list()
	for i, v in enumerate(dimtags):
		if (v[0] == gdim):
			tags.append(v[1])
	return tags	


### Model dimension
gdim = 3

flag_visu = 1

### Dimensions
corethickness = 0.019
cornerradius = 0.*0.002
corewidth = 0.065
coreheight = 0.067
coreaperturewidth = 0.015
coreapertureheight = 0.047
airgap = 0.002
corecentrallegwidth = 0.019
corelegwidth = 0.5*(corewidth-(corecentrallegwidth+2*coreaperturewidth))
airRadius = 0.075
shellRadius = 1.2*airRadius
filletRadius = 0.002

strandRadius = 0.0018
distancecenterstrands = 2*strandRadius+0.002
separationcorelayer_1 = 0.0025


### Winding parameters: double helix
R_1 = 1.2*(0.5*corecentrallegwidth+separationcorelayer_1+strandRadius) # Inner radius, first loop
R_2 = R_1+distancecenterstrands # horizontal space between strands
R_3 = strandRadius # radius of strands
a = 0.7*distancecenterstrands # vertical space between strands

nbp_1 = 5 # Number of points in inlet and outlet
nbp_2 = 100 # Number of points in helix
nbp_3 = 10 # Number of points in transition
nbt = 9 # Number of turns in helix

para = [nbp_1, nbp_2, nbp_3, nbt]

x_in = 2*R_2

# Mesh size
lc_1 = R_1 / 6

### Winding
winding = double_helix(R_1, R_2, R_3, a, x_in, para)

### Core
corebulk = occ.addBox(-0.5*corewidth, -0.5*corethickness, -0.5*coreheight, corewidth, corethickness, coreheight,  3)
coretool_1 = occ.addBox(-0.5*corecentrallegwidth-coreaperturewidth, -0.5*corethickness, -0.5*coreapertureheight, coreaperturewidth, corethickness, coreapertureheight, 4)
coretool_2 = occ.addBox(0.5*corecentrallegwidth, -0.5*corethickness, -0.5*coreapertureheight, coreaperturewidth, corethickness, coreapertureheight, 5)
coretool_3 = occ.addBox(-0.5*corewidth, -0.5*corethickness, -0.5*airgap, corewidth, corethickness, airgap,  6)
core = occ.cut([(gdim, 3)], [(gdim, i) for i in range(coretool_1,coretool_3+1)], removeObject = True, removeTool = True)
core_tags = get_gdimtags(core[0], gdim)

### Air
box_tmp = occ.addBox(-x_in, -x_in, -x_in, 2*x_in, 2*x_in, 2*x_in, 8)
air = occ.cut([(gdim, box_tmp)], winding+core[0], removeObject = True, removeTool = False)
air_tags = get_gdimtags(air[0], gdim)

### Remove duplicates
occ.removeAllDuplicates
	
### Symchronize OCCT with Gmsh
occ.synchronize()

meshing.setSize(mod.getEntities(0), lc_1)


air_surf = mod.getBoundary(air[0], combined=False, oriented=False, recursive=False)
print("air_surf", air_surf)
airsurf_tags = get_gdimtags(air_surf, gdim-1)

mod.addPhysicalGroup(gdim, [1], 1, name="Winding")
mod.addPhysicalGroup(gdim, [2, 3], 2, name="Core")
mod.addPhysicalGroup(gdim, [8], 3, name="Air")
mod.addPhysicalGroup(gdim-1, airsurf_tags[:6], 4, name="Boundary")
domain_tags = [i for i in range(2, 4)]


### Remove duplicate nodes and elements
meshing.removeDuplicateNodes
meshing.removeDuplicateElements

# ~ meshing.addHomologyRequest("Cohomology", domain_tags, [], [gdim-1])
# ~ meshing.computeHomology
	
### Colors
mod.setColor(winding, 255, 0, 0, recursive=False)
mod.setColor(core[0], 60, 60, 60, recursive=False)
mod.setColor(air[0], 50, 50, 255, recursive=False)

meshing.generate(gdim)

gmsh.write(name+".msh")
gmsh.write(name+".stp")

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

# Build the input file for the solver	
# ~ myFile0 = open(name+".par", "w")
# ~ with myFile0 as file:
	# ~ file.write("airRadius = "+str(airRadius)+";\n"+"shellRadius = "+str(shellRadius)+";\n")
	# ~ file.write("nbw = "+str(nbw)+";\n")
	# ~ file.write("nbID = "+str(k)+";\n")
	# ~ file.write("corethickness = 0.01;\n")
	# ~ file.write("coreID = 1;\n")
	# ~ file.write("airID = 2;\n")
	# ~ file.write("shellID = 3;\n")
	# ~ file.write("boundaryID = 4;\n")
	# ~ for i, v in enumerate(wiresp_tags):
		# ~ file.write("wirepID_"+str(i+1)+" = "+str(v)+";\n")
	# ~ for i, v in enumerate(wiresn_tags):
		# ~ file.write("wirenID_"+str(i+1)+" = "+str(v)+";\n")
	# ~ for i, v in enumerate(edgewiresp_tags):
		# ~ file.write("edgewirepID_"+str(i+1)+" = "+str(v)+";\n")
	# ~ for i, v in enumerate(edgewiresn_tags):
		# ~ file.write("edgewirenID_"+str(i+1)+" = "+str(v)+";\n")
# ~ myFile0.close()

