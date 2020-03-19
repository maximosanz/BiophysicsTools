import argparse
import numpy as np
import MDAnalysis as mda
import sys
import warnings
import os

def do_ang(u,v,returnCos=False,inRad=True):
	cosang = (u*v).sum(-1)/np.linalg.norm(u,axis=-1)/np.linalg.norm(v,axis=-1)
	cosang = np.clip(cosang,-1.,1.)
	if returnCos:
		return cosang
	ang = np.arccos(cosang)
	if not inRad:
		ang *= 180./np.pi
	return ang

def calc_atan2Rho(u,v,B0):
	vxu = np.cross(v,u)
	b0xu = np.cross(B0,u)
	dotrho = (vxu*b0xu).sum(-1)
	unorm = u / np.expand_dims(np.linalg.norm(u,axis=-1),-1)
	detrho = (unorm*np.cross(vxu,b0xu)).sum(-1)
	atan2rho = np.arctan2(detrho,dotrho)
	return atan2rho

def calc_CSA(XYZ,magn):

	XYZ = XYZ.T

	ENx = XYZ[0,0]
	ENy = XYZ[1,0]
	ENz = XYZ[2,0]

	Hx = XYZ[0,1]
	Hy = XYZ[1,1]
	Hz = XYZ[2,1]

	Cx = XYZ[0,2]
	Cy = XYZ[1,2]
	Cz = XYZ[2,2]

	magn = magn.T
	magn11 = magn[0]
	magn22 = magn[1]
	magn33 = magn[2]
	
	CSA = (((-0.29237173069689565*(-ENz + Hz))/np.sqrt((-ENx + Hx)**2 + (-ENy + Hy)**2 + (-ENz + Hz)**2) + 
		(0.9563047480219378*(Cz*(ENx)**2 + Cz*(ENy)**2 - Cx*ENx*ENz - Cy*ENy*ENz - 2*Cz*ENx*Hx + Cx*ENz*Hx + 
		ENx*ENz*Hx + Cz*(Hx)**2 - ENz*(Hx)**2 - 2*Cz*ENy*Hy + Cy*ENz*Hy + ENy*ENz*Hy + Cz*(Hy)**2 - 
		ENz*(Hy)**2 + Cx*ENx*Hz - (ENx)**2*Hz + Cy*ENy*Hz - (ENy)**2*Hz - Cx*Hx*Hz + ENx*Hx*Hz - 
		Cy*Hy*Hz + ENy*Hy*Hz))/ np.sqrt((Cz*(ENx)**2 + Cz*(ENy)**2 - Cx*ENx*ENz - Cy*ENy*ENz - 2*Cz*ENx*Hx + 
		Cx*ENz*Hx + ENx*ENz*Hx + Cz*(Hx)**2 - ENz*(Hx)**2 - 2*Cz*ENy*Hy + Cy*ENz*Hy + ENy*ENz*Hy + 
		Cz*(Hy)**2 - ENz*(Hy)**2 + Cx*ENx*Hz - (ENx)**2*Hz + Cy*ENy*Hz - (ENy)**2*Hz - Cx*Hx*Hz + ENx*Hx*Hz - 
		Cy*Hy*Hz + ENy*Hy*Hz)**2 + (-(Cy*ENx*ENy) + Cx*(ENy)**2 - Cz*ENx*ENz + Cx*(ENz)**2 + Cy*ENy*Hx - 
		(ENy)**2*Hx + Cz*ENz*Hx - (ENz)**2*Hx + Cy*ENx*Hy - 2*Cx*ENy*Hy + ENx*ENy*Hy - Cy*Hx*Hy + ENy*Hx*Hy + 
		Cx*(Hy)**2 - ENx*(Hy)**2 + Cz*ENx*Hz - 2*Cx*ENz*Hz + ENx*ENz*Hz - Cz*Hx*Hz + ENz*Hx*Hz + Cx*(Hz)**2 - 
		ENx*(Hz)**2)**2 + (Cy*(ENx)**2 - Cx*ENx*ENy - Cz*ENy*ENz + Cy*(ENz)**2 - 2*Cy*ENx*Hx + Cx*ENy*Hx + 
		ENx*ENy*Hx + Cy*(Hx)**2 - ENy*(Hx)**2 + Cx*ENx*Hy - (ENx)**2*Hy + Cz*ENz*Hy - (ENz)**2*Hy - Cx*Hx*Hy + 
		ENx*Hx*Hy + Cz*ENy*Hz - 2*Cy*ENz*Hz + ENy*ENz*Hz - Cz*Hy*Hz + ENz*Hy*Hz + Cy*(Hz)**2 - 
		ENy*(Hz)**2)**2))**2*magn11 + ((-(Cy*ENx) + Cx*ENy + Cy*Hx - ENy*Hx - Cx*Hy + ENx*Hy)**2*magn22) / ((-(Cy*ENx) + 
		Cx*ENy + Cy*Hx - ENy*Hx - Cx*Hy + ENx*Hy)**2 + (Cz*ENx - Cx*ENz - Cz*Hx + ENz*Hx + 
		Cx*Hz - ENx*Hz)**2 + (-(Cz*ENy) + Cy*ENz + Cz*Hy - ENz*Hy - Cy*Hz + ENy*Hz)**2) + 
		((0.9563047480219378*(-ENz + Hz))/np.sqrt((-ENx + Hx)**2 + (-ENy + Hy)**2 + (-ENz + Hz)**2) + 
		(0.29237173069689565*(Cz*(ENx)**2 + Cz*(ENy)**2 - Cx*ENx*ENz - Cy*ENy*ENz - 2*Cz*ENx*Hx + Cx*ENz*Hx + 
		ENx*ENz*Hx + Cz*(Hx)**2 - ENz*(Hx)**2 - 2*Cz*ENy*Hy + Cy*ENz*Hy + ENy*ENz*Hy + Cz*(Hy)**2 - ENz*(Hy)**2 + 
		Cx*ENx*Hz - (ENx)**2*Hz + Cy*ENy*Hz - (ENy)**2*Hz - Cx*Hx*Hz + ENx*Hx*Hz - Cy*Hy*Hz + 
		ENy*Hy*Hz)) / np.sqrt((Cz*(ENx)**2 + Cz*(ENy)**2 - Cx*ENx*ENz - Cy*ENy*ENz - 2*Cz*ENx*Hx + Cx*ENz*Hx + 
		ENx*ENz*Hx + Cz*(Hx)**2 - ENz*(Hx)**2 - 2*Cz*ENy*Hy + Cy*ENz*Hy + ENy*ENz*Hy + Cz*(Hy)**2 - ENz*(Hy)**2 + 
		Cx*ENx*Hz - (ENx)**2*Hz + Cy*ENy*Hz - (ENy)**2*Hz - Cx*Hx*Hz + ENx*Hx*Hz - Cy*Hy*Hz + ENy*Hy*Hz)**2 + 
		(-(Cy*ENx*ENy) + Cx*(ENy)**2 - Cz*ENx*ENz + Cx*(ENz)**2 + Cy*ENy*Hx - (ENy)**2*Hx + Cz*ENz*Hx - (ENz)**2*Hx + 
		Cy*ENx*Hy - 2*Cx*ENy*Hy + ENx*ENy*Hy - Cy*Hx*Hy + ENy*Hx*Hy + Cx*(Hy)**2 - ENx*(Hy)**2 + Cz*ENx*Hz - 2*Cx*ENz*Hz + 
		ENx*ENz*Hz - Cz*Hx*Hz + ENz*Hx*Hz + Cx*(Hz)**2 - ENx*(Hz)**2)**2 + (Cy*(ENx)**2 - Cx*ENx*ENy - Cz*ENy*ENz + 
		Cy*(ENz)**2 - 2*Cy*ENx*Hx + Cx*ENy*Hx + ENx*ENy*Hx + Cy*(Hx)**2 - ENy*(Hx)**2 + Cx*ENx*Hy - (ENx)**2*Hy + 
		Cz*ENz*Hy - (ENz)**2*Hy - Cx*Hx*Hy + ENx*Hx*Hy + Cz*ENy*Hz - 2*Cy*ENz*Hz + ENy*ENz*Hz - Cz*Hy*Hz + 
		ENz*Hy*Hz + Cy*(Hz)**2 - ENy*(Hz)**2)**2))**2*magn33).T
	
	return CSA


parser = argparse.ArgumentParser(description='Calculate solid-state NMR observables from molecular dynamics trajectories')

parser.add_argument('-s',metavar='structFile',help='Protein structure input file (GRO format)',required=True)
parser.add_argument('-f',metavar='trajFile',help='Trajectory input file',required=True)
parser.add_argument('-ssnmr',metavar='ssnmr',help='File with one line per observable, following the GROMACS-ssNMR format',required=True)
parser.add_argument('-multi',metavar='M',help='Run the script with mpirun in parallel for M replicas',default='0')
parser.add_argument('-kDC',metavar='K',help='Use K by default as a DC constant',default='1.0')
parser.add_argument('-periodicity',metavar='P',help='Number of atoms in each periodic unit (e.g. protomer)',default="1")
parser.add_argument('-nmolecules',metavar='M',help='Number of identical molecules (e.g. protomers)',default="1")
parser.add_argument('-o',metavar='Rho',help='File with Q-factors during the simulation',default=None)
parser.add_argument('-BC',metavar='BackCalc',help='File with back-calculated observables',default=None)

args = parser.parse_args()

script_dir = os.path.dirname(os.path.realpath(__file__))
RESTRAINT_TYPES = 4
N_GEOMETRIES = 6
DC_Geometries = np.array([ np.loadtxt(script_dir+"/GEOMETRIES/Geometry_{}.txt".format(i)) for i in range(N_GEOMETRIES) ])
N_Grid = 200

TauEdges = np.linspace(0,np.pi,N_Grid+1)
RhoEdges = np.linspace(-np.pi,np.pi,N_Grid+1)

TAU = ( TauEdges[1:] + TauEdges[:-1] ) / 2.
RHO = ( RhoEdges[1:] + RhoEdges[:-1] ) / 2.

B0 = np.array([0.,0.,1.])

multi = int(args.multi)
rank = 0
trajfn = args.f
if multi:
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	trajfn_suf = "."+trajfn.split('.')[-1]
	trajfn = trajfn[:-len(trajfn_suf)]+str(rank)+trajfn_suf

if not rank:
	o = open(args.o,'w')

u = mda.Universe(args.s,trajfn)

N_Molecules = int(args.nmolecules)
periodicity = int(args.periodicity)

input_ls = open(args.ssnmr).readlines()

N_Restr = len(input_ls)

R_Type = np.zeros(N_Restr,dtype=int)
R_Atoms = np.zeros((N_Restr,N_Molecules,3),dtype=int)
R_Exp = np.zeros(N_Restr)
R_Err = np.zeros(N_Restr)
R_kDC = np.ones(N_Restr)
R_kDC[:] = float(args.kDC)
R_CSA_magn = np.zeros((N_Restr,3))
R_CSA_magn[:] = np.array([64.0,76.0,216.9])
R_avg_eq = np.zeros(N_Restr,dtype=int)
R_geom = np.zeros(N_Restr,dtype=int)
R_Calc = np.zeros((N_Restr,N_Molecules))
R_FinalCalc = np.zeros((N_Restr))

# ssnmr file parsing
for r,l in enumerate(input_ls):
	c = l.split()
	tp = int(c[0])
	R_Type[r] = tp
	ats = np.zeros((N_Molecules,3),dtype=int)
	for i in range(3):
		if i == 2 and tp not in [0,3]:
			break
		at = int(c[i+1])
		ats[:,i] = np.arange(at,at+periodicity*N_Molecules,periodicity) - 1
		R_Atoms[r] = ats
	for idx,arr in { 3 : R_Exp, 4 : R_Err }.items():
		if tp in [0,3]:
			idx += 1
		if len(c) > idx:
			arr[r] = float(c[idx])
	if not tp and len(c) > 8:
		R_CSA_magn[r] = np.array(c[6:9],dtype=float)
	elif tp in [1,3]:
		idx = 5
		if tp == 3:
			idx += 1
		if len(c) > idx:
			R_kDC[r] = float(c[idx])
		if tp == 3 and len(c) > 7:
			R_geom[r] = float(c[7])
	avg_eq_pos = { 0 : 9 , 1 : 6 , 3 : 8 }
	idx = avg_eq_pos[tp]
	if len(c) > idx:
		try:
			R_avg_eq[r] = int(c[idx])
		except ValueError:
			pass

Traj = u.trajectory
Nframes = len(Traj)

for ts in Traj:

	if not rank:
		perc = int(ts.frame*100/Nframes)
		print("\r  >> Processed_frames: {}/{}  ({}%)    ".format(ts.frame,Nframes,perc), end='')
		sys.stdout.flush()

	Coords = ts.positions[R_Atoms]

	for R_t in range(RESTRAINT_TYPES):
		# Distances not implemented here
		if R_t == 2:
			continue
		R_IDX = R_Type == R_t
		R_Coords = Coords[R_IDX]
		if R_t == 0:
			CSA = calc_CSA(R_Coords,R_CSA_magn[R_IDX])
			R_Calc[R_IDX] = CSA
		elif R_t == 1:
			v = R_Coords[:,:,1] - R_Coords[:,:,0]
			DC = np.expand_dims(R_kDC[R_IDX],1) * 0.5 * ( 3 * do_ang(v,B0,returnCos=True)**2 - 1. )
			R_Calc[R_IDX] = DC
		elif R_t == 3:
			u = R_Coords[:,:,0] - R_Coords[:,:,1]
			v = R_Coords[:,:,2] - R_Coords[:,:,1]
			Tau = do_ang(u,B0)
			Rho = calc_atan2Rho(u,v,B0)
			TauIDX = np.round( Tau / np.pi * (N_Grid-1) ).astype(int)
			RhoIDX = np.round( (Rho + np.pi) / (2.0*np.pi) * (N_Grid-1) ).astype(int)
			GeomIDX = np.expand_dims(R_geom[R_IDX],-1)
			DC = np.expand_dims(R_kDC[R_IDX],1) * DC_Geometries[GeomIDX,TauIDX,RhoIDX]
			R_Calc[R_IDX] = DC

	R_CalcAvg = R_Calc.mean(-1)

	for i in range(np.max(R_avg_eq)+1):
		if not i:
			continue
		Avg_IDX = np.where(R_avg_eq == i)
		R_CalcAvg[Avg_IDX] = np.mean(R_CalcAvg[Avg_IDX])

	if multi:
		if not rank:
			R_RepCalc = np.zeros_like(R_CalcAvg)
		else:
			R_RepCalc = None
		comm.Reduce(
			[R_CalcAvg, MPI.DOUBLE],
			[R_RepCalc, MPI.DOUBLE],
			op = MPI.SUM,
			root = 0
		)
	else:
		R_RepCalc = R_CalcAvg
	if not rank:
		if multi:
			R_RepCalc /= float(multi)
		R_FinalCalc += R_RepCalc

		Q = np.zeros(RESTRAINT_TYPES)
		o.write("Time: {:10.3f} ; ".format(ts.time))
		for R_t in range(RESTRAINT_TYPES):
			if R_t == 2:
				continue
			R_IDX = R_Type == R_t
			Q[R_t] = np.sqrt(((R_Exp[R_IDX] - R_RepCalc[R_IDX])**2).sum() / (R_Exp[R_IDX]**2).sum())
			o.write("Q-factor {} = {:5.3f} ; ".format(R_t,Q[R_t]))
		o.write("\n")

if not rank:
	BC = open(args.BC,'w')
	R_FinalCalc /= Nframes
	for i in range(N_Restr):
		BC.write("{:6.3f}\n".format(R_FinalCalc[i]))
