import argparse
import numpy as np
import MDAnalysis as mda
import sys
import warnings

parser = argparse.ArgumentParser(description='Evaluate contacts between groups from a trajectory')

parser.add_argument('-s',metavar='structFile',help='Protein structure input file (GRO format)',required=True)
parser.add_argument('-f',metavar='trajFile',help='Trajectory input file',required=True)
parser.add_argument('-sel',metavar='SelFile',help='File with two lines - each one with an atom selection',required=True)
parser.add_argument('-epsilon',metavar='0.5',help='Threshold for contact counting (Angstrom)',default='5')
parser.add_argument('-multi',metavar='M',help='Run the script with mpirun in parallel with M processors (beware of memory usage)',default='0')
parser.add_argument('-rings',help='Compute ring-current contacts between group 1 and selected aromatic rings',default=False,action="store_true")
parser.add_argument('-residues',help='Compute the residue-wise contacts rather than atom-wise',default=False,action="store_true")
parser.add_argument('-timewise_mat',help='Add a time dimension to the matrix - be careful, this can make your memory explode!!!',default=False,action="store_true")
parser.add_argument('-collapse_mat',help='Collapse matrices onto firs axis (useful e.g. for protein-water contacts)',default=False,action="store_true")
parser.add_argument('-ring_angle',metavar='45',help='Max angle (degrees) to consider ring effect',default='45')
parser.add_argument('-cmat',metavar='CMAT.npy',help='Output contact matrix (numpy format)',required=True)
args = parser.parse_args()

multi = int(args.multi)
rank = 0
if multi:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

eps = float(args.epsilon)
ring_ang_epsilon = float(args.ring_angle)

u = mda.Universe(args.s,args.f)

sel_str = [ l.strip() for l in open(args.sel).readlines() ]

if len(sel_str) < 2:
	raise ValueError("Please include at least two lines in your -sel file")
elif len(sel_str) > 2 and not args.rings:
	warnings.warn("Option -rings not selected: Reading only the first two atom selection lines")

Traj = u.trajectory
Nframes = len(Traj)

sels = [ u.select_atoms(s) for s in sel_str ]
if not args.rings:
	sels = sels[:2]

N_ats = tuple([ len(s) for s in sels ])

if args.residues:
	Atoms = u.atoms
	Res_unique = []
	Res_counts = []
	Res_reduceIDX = []
	Res_Mapping = []
	for i,sel in enumerate(sels):
		if args.rings and i:
			break
		res = np.array([ at.resid for at in sel ])
		Res_Mapping.append(res)
		u ,c = np.unique(res,return_counts=True)
		Res_unique.append(u)
		Res_counts.append(c)
		idx = np.cumsum(c)
		idx[1:] = idx[:-1]
		idx[0] = 0
		Res_reduceIDX.append(idx)
	N_res = tuple([ len(u) for u in Res_unique ])
	
cmat_shape = N_ats
if args.residues:
	cmat_shape = N_res
if args.rings:
	N_rings = len(sels[1:])
	cmat_shape = (cmat_shape[0],N_rings)
	ring_normals = np.zeros((N_rings,3))
if args.collapse_mat:
	cmat_shape = (cmat_shape[0],)
if args.timewise_mat:
	cmat_shape = (Nframes,*cmat_shape)

cmat = np.zeros(cmat_shape,dtype=int)

cts_shape = N_ats
if args.rings:
	cts_shape = (N_ats[0],N_rings)
	ring_coseps = np.absolute(np.cos(ring_ang_epsilon/180.*np.pi))
	print("COS EPS",ring_coseps)

N_assigned = Nframes
frame0 = 0
if multi:
	N_assigned = int(Nframes/multi)
	frame0 = N_assigned*rank
	if rank == multi-1:
		N_assigned += Nframes % multi

Processed_frames = 0

for ts in Traj[frame0:frame0+N_assigned]:
	
	if multi:
		Processed_frames_mpi = comm.gather(Processed_frames, root=0)
		if not rank:
			Frame_ct = sum(Processed_frames_mpi)
	else:
		Frame_ct = Processed_frames

	if not rank:
		perc = int(Frame_ct*100/Nframes)
		print("\r  >> Processed_frames: {}/{}  ({}%)    ".format(Frame_ct,Nframes,perc), end='')
		sys.stdout.flush()
	
	coords = [ s.positions for s in sels ]

	if args.rings:
		for c_i,c in enumerate(coords[1:]):
			centroid = c.mean(0)
			coords[1+c_i] = centroid
			ring_origin = c - centroid
			uu, dd, vv = np.linalg.svd(ring_origin)
			# The third principal component is orthogonal to the ring
			ring_normals[c_i,:] = vv[2]
		coords[1] = np.array(coords[1:])

	# This is a masking approach that speeds things up if doing water contacts (many atoms on group 2) 
	lims = np.zeros((2,3))
	lims[0] = coords[0].min(0)-eps
	lims[1] = coords[0].max(0)+eps
	mask = ((coords[1] > lims[0]) & (coords[1] < lims[1])).all(-1)
	sliced_xyz = coords[1][mask]

	v = np.expand_dims(coords[0],1) - np.expand_dims(sliced_xyz,0)
	d = np.linalg.norm(v,axis=-1)
	cts_idx = np.array(d<eps)
	cts = np.array(cts_idx,dtype=int)

	if args.rings:
		# Check the angle formed by the ring normal and the relative position
		dot = (v*np.expand_dims(ring_normals,0)).sum(-1)
		cos = dot / np.linalg.norm(v,axis=-1)
		ang_cts = np.array(np.abs(cos)>ring_coseps,dtype=int)
		cts *= ang_cts

	if args.residues:
		Unmasked_cts = np.zeros(cts_shape)
		Unmasked_cts[:,mask] = cts
		for dim in range(len(cts_shape)):
			if args.rings and dim:
				break
			Unmasked_cts = np.maximum.reduceat(Unmasked_cts,Res_reduceIDX[dim],dim)
		cts = np.array(Unmasked_cts,dtype=int)
		
	if args.collapse_mat:
		cts = cts.sum(1)
	elif not args.residues:
		Unmasked_cts = np.zeros(cts_shape)
		Unmasked_cts[:,mask] = cts
		cts = np.array(Unmasked_cts,dtype=int)

	if not args.timewise_mat:
		cmat += cts
	else:
		cmat[ts.frame] += cts

	Processed_frames += 1 


if multi:
	if not rank:
		totals = np.zeros_like(cmat)
	else:
		totals = None
	comm.Reduce(
		[cmat, MPI.DOUBLE],
		[totals, MPI.DOUBLE],
		op = MPI.SUM,
		root = 0
	)
	if not rank:
		cmat = totals

if not rank:
	np.save(args.cmat,cmat)
