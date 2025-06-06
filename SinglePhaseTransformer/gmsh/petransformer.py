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
gmsh.option.setColor("General.BackgroundGradient", 255, 255, 255)
# ~ gmsh.option.setColor("General.Background", 255, 255, 255)
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
airRadius = 0.075
shellRadius = 1.2*airRadius

strandRadius = 0.0018
nbturns_1 = 9
distancecenterstrands = 2*strandRadius+0.001
layer1height = (nbturns_1-1) * distancecenterstrands
separationcorelayer_1 = 0.0025
xposstrand_1p = 0.5*corecentrallegwidth+separationcorelayer_1+strandRadius
yposstrand_1 = -0.5*layer1height

nbturns_2 = 8
separationlayer_1a2 = 0.001
xposstrand_2p = xposstrand_1p+separationlayer_1a2+2*strandRadius
layer2height = (nbturns_2-1) * distancecenterstrands

xposstrand_1n = (-1)*xposstrand_1p
xposstrand_2n = (-1)*xposstrand_2p

nbw = nbturns_1+nbturns_2

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
	### Core
	corebulk = cascade.addRectangle(-0.5*corewidth, -0.5*coreheight, 0., corewidth, coreheight, 1, cornerradius)
	coretool_1 = cascade.addRectangle(-0.5*corecentrallegwidth-coreaperturewidth, -0.5*coreapertureheight, 0., coreaperturewidth, coreapertureheight, 2, cornerradius)
	coretool_2 = cascade.addRectangle(0.5*corecentrallegwidth, -0.5*coreapertureheight, 0., coreaperturewidth, coreapertureheight, 3, cornerradius)
	coretool_3 = cascade.addRectangle(-0.5*corewidth, -0.5*airgap, 0., corewidth, airgap, 4, 0.*cornerradius)
	core = cascade.cut([(gdim, 1)], [(gdim, i) for i in range(coretool_1,coretool_3+1)], removeObject = True, removeTool = True)
	core_tags = get_gdimtags(core[0], gdim)
	
	### Turns
	strands_1p = list()
	for i in range(nbturns_1):
		strands_1p.append(cascade.addDisk(xposstrand_1p, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(core_tags)+i+1))

	strands_2p = list()
	for i in range(nbturns_2):
		strands_2p.append(cascade.addDisk(xposstrand_2p, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(strands_1p)+i+1))
	strands_p = strands_1p+strands_2p

	strands_1n = list()
	for i in range(nbturns_1):
		strands_1n.append(cascade.addDisk(xposstrand_1n, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(strands_2p)+i+1))

	strands_2n = list()
	for i in range(nbturns_2):
		strands_2n.append(cascade.addDisk(xposstrand_2n, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(strands_1n)+i+1))
	strands_n = strands_1n+strands_2n
	
	strands = strands_p+strands_n
	
	airbulk = cascade.addDisk(0., 0., 0., airRadius, airRadius)
	air = cascade.cut([(gdim, airbulk)], core[0]+[(gdim, i) for i in range(strands[0], strands[-1]+1)], removeObject = True, removeTool = False)
	air_tags = get_gdimtags(air[0], gdim)
	
	shellbulk = cascade.addDisk(0., 0., 0., shellRadius, shellRadius)
	shell = cascade.cut([(gdim, shellbulk)], [(gdim, i) for i in air_tags]+ \
		[(gdim, i) for i in core_tags]+[(gdim, i) for i in strands], \
		removeObject = True, removeTool = False)
	shell_tags = get_gdimtags(shell[0], gdim)
	
	cascade.removeAllDuplicates
	
	cascade.synchronize()
	
	meshing.setSize(mod.getEntities(0), lc_1)
	shell_surfNodes = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=True)
	meshing.setSize(shell_surfNodes, lc_2)
	
	shell_surf = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=False)
	shell_dt = get_gdimtags(shell_surf, gdim-1)
	
	strandsp_surfNodes = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_p)], combined=True, oriented=False, recursive=False)
	strandsp_surfNodes_tags = get_gdimtags(strandsp_surfNodes, gdim-1)
	strandsn_surfNodes = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_n)], combined=True, oriented=False, recursive=False)
	strandsn_surfNodes_tags = get_gdimtags(strandsn_surfNodes, gdim-1)

	meshing.removeDuplicateNodes
	meshing.removeDuplicateElements
	
	k = 0
	mod.addPhysicalGroup(gdim, core_tags, k+1, name="core")
	mod.addPhysicalGroup(gdim, air_tags, k+2, name="air")
	mod.addPhysicalGroup(gdim, [shell_tags[0]], k+3, name="shell")
	mod.addPhysicalGroup(gdim-1, [shell_dt[1]], k+4, name="boundary")
	
	domain_tags = list()
	for i in range(k+1,k+4):
		domain_tags.append(i)
	
	k = k+4
	wiresp_tags = list()
	wiresn_tags = list()
	for i, v in enumerate(strands_p):
		mod.addPhysicalGroup(gdim, [v], k+i+1, name="wirep_"+str(i+1))
		wiresp_tags.append(k+i+1)
	
	k = k+i+1
	for i, v in enumerate(strands_n):
		mod.addPhysicalGroup(gdim, [v], k+i+1, name="wiren_"+str(i+1))
		wiresn_tags.append(k+i+1)
	wires_tags = wiresp_tags+wiresn_tags
	
	k = k+i+1
	edgewiresp_tags = list()
	for i, v in enumerate(strandsp_surfNodes_tags):
		mod.addPhysicalGroup(gdim-1, [v], k+i+1, name="edgewirep_"+str(i+1))
		edgewiresp_tags.append(k+i+1)
	
	k = k+i+1
	edgewiresn_tags = list()
	for i, v in enumerate(strandsn_surfNodes_tags):
		mod.addPhysicalGroup(gdim-1, [v], k+i+1, name="edgewiren_"+str(i+1))
		edgewiresn_tags.append(k+i+1)
	k = k+i+1
	
	domain_tags = domain_tags+wires_tags
	for i, v in enumerate(wires_tags):
		idx = domain_tags.index(v)
		meshing.addHomologyRequest("Cohomology", domain_tags[:idx]+domain_tags[idx+1:], [], [gdim-1])
	meshing.computeHomology
	
	### Colors
	mod.setColor(core[0], 60, 60, 60, recursive=False)
	mod.setColor(air[0], 50, 50, 255, recursive=False)
	mod.setColor([shell[0][0]], 0, 100, 200, recursive=False)
	mod.setColor([(gdim-1, shell_dt[1])], 150, 0, 75, recursive=False)
	mod.setColor([(gdim, i) for i in range(strands_p[0], strands_p[-1]+1)], 255, 0, 0, recursive=False)
	mod.setColor([(gdim, i) for i in range(strands_n[0], strands_n[-1]+1)], 255, 150, 50, recursive=False)
	
	meshing.generate(gdim)
	
	gmsh.write(name+".msh")

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

# ~ # Build the input file for the solver	
myFile0 = open(name+".par", "w")
with myFile0 as file:
	file.write("airRadius = "+str(airRadius)+";\n"+"shellRadius = "+str(shellRadius)+";\n")
	file.write("nbw = "+str(nbw)+";\n")
	file.write("nbID = "+str(k)+";\n")
	file.write("corethickness = 0.01;\n")
	file.write("coreID = 1;\n")
	file.write("airID = 2;\n")
	file.write("shellID = 3;\n")
	file.write("boundaryID = 4;\n")
	for i, v in enumerate(wiresp_tags):
		file.write("wirepID_"+str(i+1)+" = "+str(v)+";\n")
	for i, v in enumerate(wiresn_tags):
		file.write("wirenID_"+str(i+1)+" = "+str(v)+";\n")
	for i, v in enumerate(edgewiresp_tags):
		file.write("edgewirepID_"+str(i+1)+" = "+str(v)+";\n")
	for i, v in enumerate(edgewiresn_tags):
		file.write("edgewirenID_"+str(i+1)+" = "+str(v)+";\n")
myFile0.close()

