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
opt.setNumber('Mesh.OptimizeThreshold', 0.3)

### Model dimension
gdim = 3

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
filletRadius = 0.002

strandRadius = 0.0018


### Winding parameters: double helix
R_1 = 0.1 # Inner radius, first loop
R_2 = 1.2*R_1 # horizontal space between strands
R_3 = 0.8*((R_2-R_1) / 2) # radius of strands
a = 0.25*R_1 # vertical space between strands

nbp_1 = 10 # Number of points in inlet and outlet
nbp_2 = 100 # Number of points in helix
nbp_3 = 50 # Number of points in transition
nbt = 3 # Number of turns in helix

para = [nbp_1, nbp_2, nbp_3, nbt]

# Mesh size
lc_1 = R_1 / 10

winding = double_helix(R_1, R_2, R_3, a, para)

x_in = 2*R_1
box_tmp = occ.addBox(-x_in, -x_in, -x_in, 2*x_in, 2*x_in, 2*x_in, 2)
box = occ.cut([(gdim, box_tmp)], winding, removeObject = True, removeTool = False)


	# ~ ### Air
	# ~ airbulk = cascade.addDisk(0., 0., 0., airRadius, airRadius)
	# ~ air = cascade.cut([(gdim, airbulk)], core[0]+[(gdim, i) for i in range(strands[0], strands[-1]+1)], removeObject = True, removeTool = False)
	# ~ air_tags = get_gdimtags(air[0], gdim)
	
	# ~ ### Shell
	# ~ shellbulk = cascade.addDisk(0., 0., 0., shellRadius, shellRadius)
	# ~ shell = cascade.cut([(gdim, shellbulk)], [(gdim, i) for i in air_tags]+ \
		# ~ [(gdim, i) for i in core_tags]+[(gdim, i) for i in strands], \
		# ~ removeObject = True, removeTool = False)
	# ~ shell_tags = get_gdimtags(shell[0], gdim)
	
	
	# ~ ### Remove duplicates
occ.removeAllDuplicates
	
	# ~ ### Symchronize OCCT with Gmsh
occ.synchronize()

meshing.setSize(mod.getEntities(0), lc_1)


	# ~ ### Boundaries (gdim-1)
	# ~ shell_surf = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=False)
	# ~ shell_dt = get_gdimtags(shell_surf, gdim-1)
	
	# ~ strandsp_surf = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_p)], combined=True, oriented=False, recursive=False)
	# ~ strandsp_surf_tags = get_gdimtags(strandsp_surf, gdim-1)
	# ~ strandsn_surf = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_n)], combined=True, oriented=False, recursive=False)
	# ~ strandsn_surf_tags = get_gdimtags(strandsn_surf, gdim-1)

	# ~ ### Mesh sizes (gdim-2)
	# ~ meshing.setSize(mod.getEntities(0), lc_1)
	# ~ shell_surfNodes = mod.getBoundary([shell[0][0]], combined=True, oriented=False, recursive=True)
	# ~ meshing.setSize(shell_surfNodes, lc_2)

	# ~ strandsp_surfNodes = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_p)], combined=True, oriented=False, recursive=True)
	# ~ meshing.setSize(strandsp_surfNodes, lc_3)
	# ~ strandsn_surfNodes = mod.getBoundary([(gdim, v) for i, v in enumerate(strands_n)], combined=True, oriented=False, recursive=True)
	# ~ meshing.setSize(strandsn_surfNodes, lc_3)
	
	# ~ core_surfNodes = mod.getBoundary(core[0], combined=True, oriented=False, recursive=True)
	# ~ core_edgeNodes = mod.getBoundary(core_surfNodes, combined=True, oriented=False, recursive=True)
	# ~ meshing.setSize(core_edgeNodes, lc_4)
	# ~ core_surfNodes_tags = get_gdimtags(core_surfNodes, gdim-1)
	# ~ core_edgeNodes_tags = get_gdimtags(core_edgeNodes, gdim-2)
	
	# ~ meshing.field.add("Distance", 1)
	# ~ meshing.field.setNumbers(1, "PointsList", core_edgeNodes_tags)
	# ~ meshing.field.setNumbers(1, "CurvesList", core_surfNodes_tags)
	# ~ meshing.field.setNumber(1, "Sampling", 100)
	
	# ~ meshing.field.add("Threshold", 2)
	# ~ meshing.field.setNumber(2, "InField", 1)
	# ~ meshing.field.setNumber(2, "SizeMin", lc_4)
	# ~ meshing.field.setNumber(2, "SizeMax", lc_1)
	# ~ meshing.field.setNumber(2, "DistMin", 0.001)
	# ~ meshing.field.setNumber(2, "DistMax", 0.01)


	# ~ ### Remove duplicate nodes and elements
	# ~ meshing.removeDuplicateNodes
	# ~ meshing.removeDuplicateElements

	# ~ ### Define group of geometries and ID
	# ~ k = 0
	# ~ mod.addPhysicalGroup(gdim, core_tags, k+1, name="core")
	# ~ mod.addPhysicalGroup(gdim, air_tags, k+2, name="air")
	# ~ mod.addPhysicalGroup(gdim, [shell_tags[0]], k+3, name="shell")
	# ~ mod.addPhysicalGroup(gdim-1, [shell_dt[1]], k+4, name="boundary")
	
	# ~ domain_tags = list()
	# ~ for i in range(k+1,k+4):
		# ~ domain_tags.append(i)
	
	# ~ k = k+4
	# ~ wiresp_tags = list()
	# ~ wiresn_tags = list()
	# ~ for i, v in enumerate(strands_p):
		# ~ mod.addPhysicalGroup(gdim, [v], k+i+1, name="wirep_"+str(i+1))
		# ~ wiresp_tags.append(k+i+1)
	
	# ~ k = k+i+1
	# ~ for i, v in enumerate(strands_n):
		# ~ mod.addPhysicalGroup(gdim, [v], k+i+1, name="wiren_"+str(i+1))
		# ~ wiresn_tags.append(k+i+1)
	# ~ wires_tags = wiresp_tags+wiresn_tags
	
	# ~ k = k+i+1
	# ~ edgewiresp_tags = list()
	# ~ for i, v in enumerate(strandsp_surf_tags):
		# ~ mod.addPhysicalGroup(gdim-1, [v], k+i+1, name="edgewirep_"+str(i+1))
		# ~ edgewiresp_tags.append(k+i+1)
	
	# ~ k = k+i+1
	# ~ edgewiresn_tags = list()
	# ~ for i, v in enumerate(strandsn_surf_tags):
		# ~ mod.addPhysicalGroup(gdim-1, [v], k+i+1, name="edgewiren_"+str(i+1))
		# ~ edgewiresn_tags.append(k+i+1)
	# ~ k = k+i+1
	
	# ~ domain_tags = domain_tags+wires_tags
	# ~ for i, v in enumerate(wires_tags):
		# ~ idx = domain_tags.index(v)
		# ~ meshing.addHomologyRequest("Cohomology", domain_tags[:idx]+domain_tags[idx+1:], [], [gdim-1])
	# ~ meshing.computeHomology
	
	# ~ ### Colors
	# ~ mod.setColor(core[0], 60, 60, 60, recursive=False)
	# ~ mod.setColor(air[0], 50, 50, 255, recursive=False)
	# ~ mod.setColor([shell[0][0]], 0, 100, 200, recursive=False)
	# ~ mod.setColor([(gdim-1, shell_dt[1])], 150, 0, 75, recursive=False)
	# ~ mod.setColor([(gdim, i) for i in range(strands_p[0], strands_p[-1]+1)], 255, 0, 0, recursive=False)
	# ~ mod.setColor([(gdim, i) for i in range(strands_n[0], strands_n[-1]+1)], 255, 150, 50, recursive=False)
	
meshing.generate(gdim)

gmsh.write(name+".msh")

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

