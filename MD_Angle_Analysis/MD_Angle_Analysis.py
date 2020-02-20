import scipy
import numpy as np
import sys
import argparse
import subprocess
import re
import commands


parser = argparse.ArgumentParser(description='Calculate tilt and azimuthal angles from MD trajectory')

parser.add_argument('-s',metavar='structFile',help='Protein structure input file',required=True)
parser.add_argument('-f',metavar='trajFile',help='Trajectory input file',required=True)
parser.add_argument('-res',metavar='helixStart,helixEnd',help='Comma separated residue range of the helix')
parser.add_argument('-res2',metavar='helix2Start,helix2End',help='Comma separated residue range of a second helix, computes the inter-helical angle INSTEAD of the tilt')
parser.add_argument('-atoms',metavar='type',help='Atoms to use for the least-squares fitting',choices=['CA','backbone'],default='CA',type=str)
parser.add_argument('-globalang',metavar='global_outputFile',help='File where the global angle will be written',default=None,type=str)
parser.add_argument('-relative',metavar='relative_outputFile',help='File where the relative angles will be written',default=None,type=str)
parser.add_argument('-azimuthalF',metavar='azimuthal_outputFile',help='File where the azimuthal angles will be written',default=None,type=str)
parser.add_argument('-azimuthalRes',metavar='AZ',help='Residue that defines the azimuthal angle',default=None,type=int)
parser.add_argument('-nprotomers',metavar='N',help='Number of protomers',default='1',type=int)
parser.add_argument('-invert',help='Default vector will go from lower to higher residue number',action='store_true')
parser.add_argument('-invertB0az',help='Invert the B0 for the azimuthal calculation (yields the 180-flipped result)',action='store_true')
parser.add_argument('-idealAz',help='Calculate the azimuthal based on the whole helix (assuming ideal geometry)',action='store_true')
parser.add_argument('-o',metavar='outputFile',help='Output file',required=True)
args = parser.parse_args()

isTraj = False
isStruct = False

if args.azimuthalF != args.azimuthalRes and None in (args.azimuthalF,args.azimuthalRes):
    parser.error("Please provide both -azimuthalF and -azimuthalRes for azimuthal calculation")

if args.idealAz and args.atoms != "CA":
    parser.error("Sorry, the ideal azimuthal calculation works only with the CA option for atoms")

if args.azimuthalF is not None:
	azo = open(args.azimuthalF,'w')

if args.res is not None:
	hrange = map(int,args.res.split(','))
if args.res2 is not None:
	h2range = map(int,args.res2.split(','))

atoms = ['CA']

if args.atoms == 'backbone':
	atoms.extend(['N','C','O'])

o = open(args.o,'w')
if args.globalang is not None:
	globalo = open(args.globalang,'w')
if args.relative is not None:
	relo = open(args.relative,'w')

B0 = np.zeros((args.nprotomers,3),dtype=float)
B0[:,2] = 1.
if args.invert:
	B0[:,2] = -1.

f = open(args.s)
ls = f.readlines()[2:-1]

def get_atnos(ls,hrange,atoms=atoms,azRes=args.azimuthalRes,nprotomers=args.nprotomers):
	atnos = []
	protatnos = []
	atind = 0
	prot = 0
	if azRes is not None:
		azprotatnos = []
		azatnos = []

	for l in ls:
		resno = int(l[:5].strip())
		at = l[13:15].strip()
		if at in atoms and resno >= hrange[0] and resno <= hrange[1]:
			atno = int(l[15:20].strip())
			protatnos.append(atno)
			if azRes is not None and resno == azRes:
				azprotatnos.append(atno)
			maxat = atno

			if resno == hrange[1]:
				atind += 1
				if atind == len(atoms):
					atind = 0
					atnos.append(protatnos)
					protatnos = []
					if azRes is not None:
						azatnos.append(azprotatnos)
						azprotatnos = []
					
					prot += 1
					if prot == nprotomers:
						break
	atnos = np.array(atnos)
	if azRes is not None:
		azatnos = np.array(azatnos)
	return atnos,azatnos,maxat


atnos, azatnos, maxat = get_atnos(ls,hrange)
coords = np.zeros((atnos.shape[0],atnos.shape[1],3))

if args.res2 is not None:
	atnos2, azatnos2, maxat2 = get_atnos(ls,h2range,azRes=False)
	coords2 = np.zeros((atnos2.shape[0],atnos2.shape[1],3))


if args.azimuthalRes is not None:
	azcoords = np.zeros((azatnos.shape[0],azatnos.shape[1],3))

curt=0
trj = 0

def calc_least_SQ(coords):
	# Calculate the centroid of the points
	datamean = coords.mean(axis=1)
	# Do an SVD on the mean-centered data.
	uu, dd, vv = np.linalg.svd(coords - np.expand_dims(datamean,axis=1))
	# Correct the direction of the first PC
	direction = coords[:,-1,:]-coords[:,0,:]
	direction /= np.linalg.norm(direction,axis=-1)
	# Direction should be the same as the first PC, otherwise correct it
	reverse = np.expand_dims(np.array(np.linalg.norm(direction-vv[:,0],axis=-1) <= 1,dtype=float) * 2 - 1,-1)
	vv[:,0] *= reverse
	return vv[:,0],dd[:,0]

def calc_Az(coords,azcoords,OLS,dd,B0,invertB0az=args.invertB0az,idealAz=args.idealAz):
	datamean = coords.mean(1)
	lowV = datamean-OLS*dd/2
	if idealAz:
		shiftAz = coords-lowV
	else:
		azcentroid = azcoords.mean(1)
		shiftAz = azcentroid-lowV
		shiftAz = np.expand_dims(shiftAz,1)
	# Now project shiftAz onto helix axis
	OLS = np.expand_dims(OLS,1)
	B0 = np.expand_dims(B0,1)
	proj = np.expand_dims((shiftAz*OLS).sum(-1),-1)*OLS
	rej = shiftAz-proj
	projB0 = np.expand_dims((B0*OLS).sum(-1),-1)*OLS
	rejB0 = B0-projB0
	if invertB0az:
		rejB0 = -rejB0
	HxB0 = np.cross(OLS,B0)
	Hxrej = np.cross(OLS,rej)
	dot = (rej*rejB0).sum(-1)
	det = (OLS*np.cross(rej,rejB0)).sum(-1)
	az = np.arctan2(det,dot)/np.pi*180
	return az

# This is very slow!!
#proc = subprocess.Popen(["gmxdump",'-f',args.f],stdout=subprocess.PIPE)
#for line in iter(proc.stdout.readline,''):

# This only works in unix and is very memory intensive!!
print "Reading..."
proc = commands.getoutput("gmxdump -f "+args.f)
for line in proc.split('\n'):
	l = line.strip()
	if l[:6] == "natoms":
		t = float(l.split("time=")[1][:13])
		if t<curt:
			trj+=1
		curt=t
		o.write(str(trj)+'\t'+str(t))
		if args.azimuthalRes is not None:
			azo.write(str(trj)+'\t'+str(t))
		print "\033[K","TRAJ",trj,"TIME:",t,"\r",
		sys.stdout.flush()
	elif l[:2] == 'x[':
		atno = int(l[2:].split(']')[0])+1 # Since gmxdump is zero-indexed
		index = np.where(atnos==atno)
		if args.azimuthalRes is not None:
			azindex = np.where(azatnos==atno)
		if args.res2:
			index2 = np.where(atnos2==atno)
			if len(index2[0]):
				index2 = tuple([index2[0][0],index2[1][0]])
			else:
				index2 = []
			if index2:
				xyz = map(float,l.split('{')[1].split('}')[0].split(','))
				coords2[index2[0],index2[1],:] = xyz
				if atno == maxat2:
					B0 = calc_least_SQ(coords2)

		if len(index[0]):
			index = tuple([index[0][0],index[1][0]])
		else:
			index = []
		if args.azimuthalRes is not None:
			if len(azindex[0]):
				azindex = tuple([azindex[0][0],azindex[1][0]])
			else:
				azindex = []
		if index or (args.res2 is not None and index2):
			if index:
				xyz = map(float,l.split('{')[1].split('}')[0].split(','))
				coords[index[0],index[1],:] = xyz
			if args.azimuthalRes is not None and azindex:
				xyz = map(float,l.split('{')[1].split('}')[0].split(','))
				azcoords[azindex[0],azindex[1],:] = xyz
			if atno >= maxat and (args.res2 is None or atno >= maxat2):
				OLS,dd = calc_least_SQ(coords)
				ang = np.arccos((B0*OLS).sum(-1)/ (np.linalg.norm(B0,axis=-1) * np.linalg.norm(OLS,axis=-1)))/np.pi*180
				for i in range(ang.shape[0]):
					o.write("\t"+str(ang[i]))

				az = calc_Az(coords,azcoords,OLS,dd,B0)

				if args.idealAz:
					ideal_Offset = 100.
					NRes = az.shape[1]
					azID = args.azimuthalRes - hrange[0]
					offsets = np.arange(-azID*ideal_Offset,NRes*ideal_Offset-azID*ideal_Offset,ideal_Offset)
					az = az + offsets
					periods = np.array(az/360.,dtype=int)
					az -= 360.*periods
					az[az<0] += 360.
					azstd = az.std(-1)
					az = az.mean(-1)
				else:
					az = az[:,0]

				o.write('\n')
				if args.azimuthalRes is not None:
					for i in range(az.shape[0]):
						azo.write("\t"+str(az[i]))
					azo.write('\n')
				#globalv = np.sum(allpcs,axis=0)
				globalv = np.sum(OLS,axis=0)
				globalB0 = np.sum(B0,axis=0)
				if args.globalang is not None:
					ang = np.arccos(np.dot(globalB0, globalv) / (np.linalg.norm(globalB0) * np.linalg.norm(globalv)))/np.pi*180
					globalo.write(str(trj)+'\t'+str(t)+'\t'+str(ang)+'\n')
				if args.relative is not None:
					i=0
					relo.write(str(trj)+'\t'+str(t))
					while i < args.nprotomers:
						vv = allpcs[i,:]
						ang = np.arccos(np.dot(globalv,vv) / (np.linalg.norm(globalv) * np.linalg.norm(vv)))/np.pi*180
						relo.write("\t"+str(ang))
						i += 1
					relo.write('\n')


