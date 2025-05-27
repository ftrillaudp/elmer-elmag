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
	strands_1p = list()
	for i in range(nbturns_1):
		strands_1p.append(cascade.addDisk(xposstrand_1p, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(core_tags)+i+1))
	print("strand number", strands_1p)
	strands_2p = list()
	for i in range(nbturns_2):
		strands_2p.append(cascade.addDisk(xposstrand_2p, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(strands_1p)+i+1))
	print("strand number", strands_2p)
	strands_p = strands_1p+strands_2p
	strands_1n = list()
	for i in range(nbturns_1):
		strands_1n.append(cascade.addDisk(xposstrand_1n, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(strands_2p)+i+1))
	print("strand number", strands_1n)
	strands_2n = list()
	for i in range(nbturns_2):
		strands_2n.append(cascade.addDisk(xposstrand_2n, yposstrand_1+i*distancecenterstrands, 0., \
		strandRadius, strandRadius, max(strands_1n)+i+1))
	print("strand number", strands_2n)
	strands_n = strands_1n+strands_2n
	strands = strands_p+strands_n
	print("strands", strands)
	
	print("core[0]", core[0]+[(gdim, i) for i in range(strands[0], strands[-1]+1)])
	
	airbulk = cascade.addDisk(0., 0., 0., airRadius, airRadius)
	air = cascade.cut([(gdim, airbulk)], core[0]+[(gdim, i) for i in range(strands[0], strands[-1]+1)], removeObject = True, removeTool = False)
	air_tags = get_gdimtags(air[0], gdim)
	print("air", air)
	
	shellbulk = cascade.addDisk(0., 0., 0., shellRadius, shellRadius)
	shell = cascade.cut([(gdim, shellbulk)], [(gdim, i) for i in air_tags], removeObject = True, removeTool = False)
	shell_tags = get_gdimtags(shell[0], gdim)
	print("shell", shell_tags)
	
	cascade.removeAllDuplicates
	
	cascade.synchronize()
	
	meshing.setSize(mod.getEntities(0), lc_1)
	shell_surfNodes = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=True)
	meshing.setSize(shell_surfNodes, lc_2)
	
	shell_surf = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=False)
	shell_dt = get_gdimtags(shell_surf, gdim-1)
	
	domainnocoils_tags = list()
	coilsp_tags = list()
	coilsn_tags = list() 
	for i, v in enumerate(strands_p):
		# ~ mod.addPhysicalGroup(gdim, [v], i+1, name="coilsp_"+str(i+1))
		# ~ domain_tags.append(i+1)
		coilsp_tags.append(v)
	k = i+1
	for i, v in enumerate(strands_n):
		# ~ mod.addPhysicalGroup(gdim, [v], k+i+1, name="coilsn_"+str(i+1))
		# ~ domain_tags.append(k+i+1)
		coilsn_tags.append(v)
	k = k+i+1
	mod.addPhysicalGroup(gdim, coilsp_tags, k+1, name="coilsp")
	mod.addPhysicalGroup(gdim, coilsn_tags, k+2, name="coilsn")
	k = k+2
	mod.addPhysicalGroup(gdim, core_tags, k+1, name="core")
	domainnocoils_tags.append(k+1)
	mod.addPhysicalGroup(gdim, air_tags, k+2, name="air")
	domainnocoils_tags.append(k+2)
	mod.addPhysicalGroup(gdim, [shell_tags[0]], k+3, name="shell")
	domainnocoils_tags.append(k+3)
	mod.addPhysicalGroup(gdim-1, [shell_dt[1]], k+4, name="boundary")
	domainnocoils_tags.append(k+4)

	
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
	
	meshing.addHomologyRequest("Homology", domainnocoils_tags, coilsn_tags+coilsp_tags)

	
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

