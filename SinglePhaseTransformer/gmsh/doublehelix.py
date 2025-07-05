import gmsh
import numpy as np

mod = gmsh.model
occ = mod.occ

def double_helix(R_1, R_2, R_3, a, x_in, para):
	
	if (R_3 == ((R_2-R_1) / 2)):
		print("The radius of the strand is changed autamically")
		R_3 = 0.9*((R_2-R_1) / 2)
	
	nbp_1 = para[0]
	nbp_2 = para[1]
	nbp_3 = para[2]
	nbt_1 = para[3]
	nbt_2 = nbt_1+0.5 # Number of turns in helix 2 (outer)
	nbt_3 = 1 # One turn in the transition by default
	p = list()

	# Start from a line, build the helix and finish with a parallel line!
	# Inlet
	for i in range(nbp_1):
		occ.addPoint(R_1, ((i-nbp_1) / nbp_1) * x_in, 0., 1, 1000+i)
		p.append(1000+i)
	del p[-1]
	# Inner helix
	for i in range(nbp_2+1):
		t = i * (2. * np.pi) * (nbt_1 / nbp_2)
		occ.addPoint(R_1*np.cos(t), R_1*np.sin(t), a* i * (nbt_1 / nbp_2), 1, 2000 + i)
		p.append(2000+i)
	del p[-1]
	# Transition
	DR = R_2-R_1
	for i in range(nbp_3+1):
		t = i * (2. * np.pi) * (nbt_3 / nbp_3)
		R = (DR / np.exp(2*np.pi))*np.exp(2*np.pi * (i+1)/(nbp_3+1)) + R_1
		occ.addPoint(R*np.cos(t), R*np.sin(t), a * nbt_1 + a* i * (nbt_3 / nbp_3), 1, 3000 + i)
		p.append(3000+i)
	del p[-1]
	# Outer helix
	for i in range(nbp_2+1):
		t = i * (2. * np.pi) * (nbt_2 / nbp_2)
		occ.addPoint(R_2*np.cos(t), R_2*np.sin(t), a * (nbt_1+nbt_3) - a* i * (nbt_2 / nbp_2), 1, 4000 + i)
		p.append(4000+i)
	# Outlet
	for i in range(nbp_1+1):
		occ.addPoint(-R_2, -(i / nbp_1) * x_in, a * (nbt_1+nbt_3-nbt_2), 1, 5000+i)
		p.append(5000+i)
	del p[-nbp_1-1]
	
	occ.addBSpline(p, 1000)
	# A wire is like a curve loop, but open:
	occ.addWire([1000], 1000)

	# ~ # We define the shape we would like to extrude along the spline (a disk):
	occ.addDisk(R_1, -x_in, 0., R_3, R_3, 6000)
	occ.rotate([(2, 6000)], R_1, -x_in, 0., 1, 0, 0, np.pi / 2)
	winding = occ.addPipe([(2, 6000)], 1000, 'DiscreteTrihedron')
	occ.translate(winding, 0, 0, -0.5*(nbt_1+1)*a)
	
	# A bit of cleaning 
	occ.remove([(2, 6000), (1, 1000), (1, 1001)]+[(0, v) for i, v in enumerate(p)])
	
	return [winding, p]
