import numpy as np
import argparse
import string
from Bio.PDB import *

argparser = argparse.ArgumentParser(description='Generate an elongated amyloid fibre structure')

argparser.add_argument('-pdb',metavar='PDB_file',help='Protein structure input file (PDB format)',required=True)
argparser.add_argument('-unit',metavar='A,B',help='Comma-separated string of chains that compose a single repeating unit (use whole structure if not provided)',required=False,default=None)
argparser.add_argument('-repeated',metavar='A,C',help='Comma-separated string of two repeated chains to calculate fibre axis (Z-axis used if not provided)',required=False,default=None)
argparser.add_argument('-residue_range',metavar='1,10',help='Selection range to limit the residues used in the axis calculation',required=False,default=None)
argparser.add_argument('-n_units',metavar='10',help='Number of repeated units to include in the fibre',required=False,default='10')
argparser.add_argument('-beta_spacing',metavar='4.8',help='Inter-strand amyloid spacing (4.75 Angstrom by default)',required=False,default='4.75')
argparser.add_argument('-pitch',metavar='100',help='Helical pitch of one full rotation (nm)',default='0')
argparser.add_argument('-right_handed',help='Build a right-handed filament (they are left-handed by default)',action='store_true',default=False)
argparser.add_argument('-curvature',metavar='10',help='Include a curvature of this radius (nm)',default='0')
argparser.add_argument('-curve_axis',metavar='1.0,0.0,0.0',help='(Optional) Comma-separated vector specifying the axis to apply the curvature',default=None)
argparser.add_argument('-invert_curvature',help='Invert the direction of the curvature',action='store_true',default=False)
argparser.add_argument('-offset',metavar='N',help='Include an offset of N units (can be float)',default=0)
argparser.add_argument('-seed',metavar='S',help='Random seed',default=0)
argparser.add_argument('-o',metavar='PDB_output',help='Output PDB file',required=True)
args = argparser.parse_args()

backbone = ["CA","C","N"]
unit_chains = None
if args.unit is not None:
	unit_chains = args.unit.split(',')
N_Units = int(args.n_units)
beta_spacing = float(args.beta_spacing)
R_curvature = float(args.curvature) * 10.0
curve_axis = args.curve_axis
offset = float(args.offset)
seed = int(args.seed)
if curve_axis is not None:
	curve_axis = np.array(curve_axis.split(','),dtype=float)

residue_range = args.residue_range
if residue_range is not None:
	residue_range = tuple([ int(r) for r in args.residue_range.split(',') ])

parser = PDBParser(QUIET=True)
structure = parser.get_structure('raw', args.pdb)
model = next(structure.get_models())

pitch = float(args.pitch) * 10.0
u_rot = 0.0
if pitch:
	u_rot = 2*np.pi / (pitch / beta_spacing)
	if args.right_handed:
		u_rot = -u_rot

def rotation_matrix(axis, theta):
	"""
	Return the rotation matrix associated with counterclockwise rotation about
	the given axis by theta radians.
	"""
	axis = np.asarray(axis)
	axis = axis / np.sqrt(np.dot(axis, axis))
	a = np.cos(theta / 2.0)
	b, c, d = -axis * np.sin(theta / 2.0)
	aa, bb, cc, dd = a * a, b * b, c * c, d * d
	bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
	return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
					 [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
					 [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def get_unpacked_list(self):
	"""
	Returns all atoms from the residue,
	in case of disordered, keep only first alt loc and remove the alt-loc tag
	"""
	atom_list = self.get_list()
	undisordered_atom_list = []
	for atom in atom_list:
		if atom.is_disordered():
			atom.altloc=" "
			undisordered_atom_list.append(atom)
		else:
			undisordered_atom_list.append(atom)
	return undisordered_atom_list
Residue.Residue.get_unpacked_list = get_unpacked_list

Unit = Structure.Structure('unit')
u_model = Model.Model(0)
Unit.add(u_model)
for chain in structure.get_chains():
	if unit_chains is None or chain.id in unit_chains:
		u_model.add(chain.copy())

COM = np.array([ atom.get_coord() for atom in Unit.get_atoms() if atom.name in backbone ]).mean(0)
Unit_XYZ = np.array([ atom.get_coord() for atom in Unit.get_atoms() if atom.name in backbone ]) - COM


F_Axis = np.array([0., 0., 1.])

def pick_coords(chain,residue_range):
	xyz = []
	for r in chain.get_residues():
		resi = r.id[1]
		if residue_range is None or (resi >= residue_range[0] and resi < residue_range[1]):
			for atom in r:
				if atom.name in backbone:
					xyz.append(atom.get_coord())
	return np.array(xyz)

if args.repeated is not None:
	rep_chains = args.repeated.split(',')
	if len(rep_chains) != 2:
		raise ValueError("Only two chains can be specified to define the fibre axis")
	xyz0 = None
	xyz1 = None
	for chain in structure.get_chains():
		if chain.id == rep_chains[0]:
			xyz0 = pick_coords(chain,residue_range)
		elif chain.id == rep_chains[1]:
			xyz1 = pick_coords(chain,residue_range)
			#xyz1 = np.array([ atom.get_coord() for atom in chain.get_atoms() if atom.name in backbone ])
	if xyz0 is None or xyz1 is None or xyz0.shape != xyz1.shape:
		raise ValueError("Selected repeated chains are not identical or not found in the structure")
	F_Axis = (xyz1 - xyz0).mean(0)
	F_Axis /= np.linalg.norm(F_Axis)

if curve_axis is not None:
	displace_axis = np.cross(curve_axis,F_Axis)

if R_curvature and curve_axis is None:
	if seed:
		np.random.seed(seed)
	rd = np.random.uniform(-1,1,size=(3,))
	x = np.cross(F_Axis,rd)
	x /= np.linalg.norm(x)
	y = np.cross(F_Axis,x)
	y /= np.linalg.norm(y)
	coord2d = np.array([x,y,F_Axis]).T
	inv_coord2d = np.linalg.inv(coord2d)
	proj = Unit_XYZ @ coord2d
	proj2d = proj[:,:2]
	# Find the principal axis of variance along the cross-section of the fibre
	# This is the default axis used to curve the fibre if none is defined
	u, s, vh = np.linalg.svd(proj2d)
	curve2d = np.concatenate([vh[0],[0.]])
	displace2d = np.concatenate([vh[1],[0.]])
	# Transform back to XYZ coordinate space
	curve_axis = curve2d @ inv_coord2d
	displace_axis = displace2d @ inv_coord2d
	'''
	import matplotlib.pyplot as plt
	plt.scatter(proj2d[:,0],proj2d[:,1])
	plt.plot([-10*curve2d[0],10*curve2d[0]],[-10*curve2d[1],10*curve2d[1]],color='red')
	plt.show()
	'''

Fibre = Structure.Structure('fibre')
f_model = Model.Model(0)
Fibre.add(f_model)

XYZ = np.array([ atom.get_coord() for atom in Unit.get_atoms() ])

cs = list(string.ascii_uppercase)
cID_l = [ c for c in cs ]
cID_l += [ c1+c for c1 in cs for c in cs ]
cID = ( c for c in cID_l )

for U in range(N_Units):
	U_off = float(U) + offset
	for chain in Unit.get_chains():
		f_chain = chain.copy()
		f_chain.id = next(cID)
		for atom in f_chain.get_atoms():
			xyz = atom.get_coord()
			xyz -= COM
			U_F_Axis = np.copy(F_Axis)
			if R_curvature:
				theta_curve = (beta_spacing / R_curvature) * U_off
				if args.invert_curvature:
					theta_curve = -theta_curve
				curve_rotmat = rotation_matrix(curve_axis,theta_curve)
				U_F_Axis = F_Axis @ curve_rotmat
			if u_rot:
				rotmat = rotation_matrix(U_F_Axis, u_rot*U_off)
				xyz = xyz @ rotmat
			if R_curvature:
				xyz += R_curvature * displace_axis
				xyz = xyz @ curve_rotmat
				xyz -= R_curvature * displace_axis
			else:
				xyz += beta_spacing*U_off*F_Axis
			atom.set_coord(xyz)
		f_model.add(f_chain)

# These functions are re-defined to cope with formatting issues
# 1. When there are more than 26 chains, we start repeating labels
# 2. When there are more than 100000 atoms, the spacing is adjusted and the atom number loses meaning

def get_id0(self):
	return self.id[-1]

def get_serial_number5(self):
	return int(str(self.serial_number)[:5])

Chain.Chain.get_id = get_id0
Atom.Atom.get_serial_number = get_serial_number5

io=PDBIO()
io.set_structure(Fibre)
io.save(args.o,preserve_atom_numbering=True)
