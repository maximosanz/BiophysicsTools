import numpy as np
from decimal import Decimal

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


def do_ang(u,v,returnCos=False,inRad=True):
	cosang = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v)
	cosang = np.clip(cosang,-1.,1.)
	if returnCos:
		return cosang
	ang = np.arccos(cosang)
	if not inRad:
		ang *= 180./np.pi
	return ang

def do_sine(u,v):
	sineang = np.linalg.norm(np.cross(u,v))/np.linalg.norm(u)/np.linalg.norm(v)
	return sineang

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def midpoints(array):
    return ( array[1:] + array[:-1] ) / 2.

# Class for the different types of geometries modelled (e.g. CH2-sp2 , CH-sp3 , N-CH3x3 (choline group))

class Geometry:

	'''
	# Class constructor builds the geometry from an example set of atomic coordinates for that group and an array specifying which atoms are bonded
	# Stores all the bonds in the geometry
	# xyz is an Nx3 np.array where N is the number of atoms in the geometry
	# bonded is an NxN np.array which specifies which pairs of atoms are bonded (in the lower triangular half of the array)
	# The directionality of a bond can be reversed by specifying -1 instead of 1. By default it will go toward the atom with the higher index in the array.
	# DCbond_IDs is a list of integers containing the indices of bonds involved in the calculation of the DC
	# Tilt ID is the bond used to define the tilt angle, and azID to define the azimuthal
	# secondOrderTau is a list with bond IDs, it gives the option to compute the DCs through an externally provided TAUproj array, using the tau of those bonds in the geometry
	#
	'''
	def __init__(self, name, xyz, bonded, DCbond_IDs, tiltID, azID, 
		kDC=1.0, 
		B0=np.array([0.,0.,1.]), 
		gridsize=(200,200), 
		tiltRotv=np.array([np.sqrt(2.0),-np.sqrt(2.0),0.]),
		secondOrderTau=None):

		self.name = name
		self.xyz = xyz
		self.bonded = bonded
		self.bondsFromXyz()
		self.Nbonds = self.bonds.shape[0]
		self.DCbond_IDs = DCbond_IDs
		self.tiltID = tiltID
		self.azID = azID
		self.kDC = kDC
		self.B0 = B0
		if self.B0 is None:
			self.B0 = np.random.normal(0.,1.,size=3)
		self.gridsize = gridsize
		self.tiltRotv = tiltRotv
		# If not tiltRotv vector provided, come up with a random one orthogonal to B0
		if self.tiltRotv is None:
			rd = self.B0
			while np.allclose(rd,self.B0,atol=1e-2):
				rd = self.B0 + np.random.normal(0.,1.,size=3)
			self.tiltRotv = np.cross(self.B0,rd)
		self.costau_edges = np.linspace(-1.,1.,self.gridsize[0]+1)
		self.costau = midpoints(self.costau_edges)
		self.tau = np.arccos(self.costau)
		self.cosrho_edges = np.linspace(-1.,1.,self.gridsize[1]+1)
		self.cosrho = midpoints(self.cosrho_edges)
		self.rho = np.arccos(self.cosrho)
		self.atan2rho_edges = np.linspace(-np.pi,np.pi,self.gridsize[1]+1)
		self.atan2rho = midpoints(self.atan2rho_edges)
		self.tauANG_edges = np.linspace(0.,np.pi,self.gridsize[0]+1)
		self.tauANG = midpoints(self.tauANG_edges)
		self.DC_surface = np.zeros(self.gridsize)
		self.BC_ders = np.zeros((12,*self.gridsize))
		self.TAUproj = np.zeros(self.gridsize[0])
		self.RHOproj = np.zeros(self.gridsize[1])
		self.tilt2ID = secondOrderTau
		if self.tilt2ID is not None:
			self.TAU2proj = np.zeros(self.gridsize[0])
		return None


	def bondsFromXyz(self):
		N = self.xyz.shape[0]
		bonds = []
		for i in range(N):
			for j in range(i):
				isbonded = self.bonded[i,j]
				if not isbonded:
					continue
				bonds.append((self.xyz[i]-self.xyz[j])*isbonded)
		self.bonds = np.array(bonds)
		return self


	'''
	# Build the geometry from a Gro file
	# Must provide the filename as first argument
	# atnames_list is a list of lists, with each internal list containing a number of atom names from the GRO file.
	# The length of all internal lists must be the same, and it must match the first dimension of the "xyz" array provided in the class constructor
	# The order of atoms must also match the one used in the "bonded" array provided in the class constructor
	# IMPORTANT: The first atom of each internal list MUST be the first atom of that list to appear int the GRO file!
	'''

	def buildFromGro(self,fname,atnames_list,checkChirality=None):
		f = open(fname)
		ls = f.readlines()[2:-1]
		xyz = []
		d = {}
		Natoms = len(atnames_list[0])
		for l in ls:
			atname = l[9:15].strip()
			for atnID, atnames in enumerate(atnames_list):
				if atname in atnames:
					atID = atnames.index(atname)
					if not atID:
						xyz.append(np.zeros((Natoms,3)))
						d[atnID] = len(xyz)-1
					coords = l[20:44].strip().split()
					coords = list(map(float,coords))
					xyz[d[atnID]][atID] = coords
		xyz = np.array(xyz)
		Ngeometries = xyz.shape[0]
		Natoms = xyz.shape[1]
		allbonds = []
		cur=-1
		for i in range(Ngeometries):
			progress = int(i*100/Ngeometries)
			if not progress % 5 and progress != cur:
				print("\rBuilding the geometry",self.name,"from GRO file",fname,"...",int(progress),"%             ", end='')
				cur = progress
			self.xyz = xyz[i]
			self.bondsFromXyz()
			self.orient(self.tiltID,self.tiltRotv)
			self.reorientAz()
			self.orient(self.tiltID,self.B0)
			if np.isnan(self.bonds).any():
				continue
			# The following checks if the two bonds in checkChirality (a tuple) are inverted with respect to the first molecule visited. If so, exchange their coords
			if checkChirality is not None and allbonds:
				chirID1 = checkChirality[0]
				chirID2 = checkChirality[1]
				chirV1 = np.copy(self.bonds[chirID1])
				chirV2 = np.copy(self.bonds[chirID2])
				origV1 = allbonds[0][chirID1]
				origV2 = allbonds[0][chirID2]
				norm11 = np.linalg.norm(chirV1-origV1)
				norm12 = np.linalg.norm(chirV1-origV2)
				norm21 = np.linalg.norm(chirV2-origV1)
				norm22 = np.linalg.norm(chirV2-origV2)
				if norm12 < norm11 and norm21 < norm22:
					self.bonds[chirID2] = chirV1
					self.bonds[chirID1] = chirV2
			allbonds.append(self.bonds)
		allbonds = np.array(allbonds)
		norms = np.linalg.norm(allbonds,axis=2)
		norms = np.expand_dims(norms,axis=2)
		allbonds /= norms
		meanv = np.mean(allbonds,axis=0)
		# This tells you the divergence between vectors in the ensemble, a sort of order parameter
		print("\n/*** Constructed the geometry",self.name,"from the GRO file",fname,"***/\n")
		print("Pseudo order parameter of bonds in the GRO file:\n",np.linalg.norm(meanv,axis=1),'\n')
		meanv /= np.expand_dims(np.linalg.norm(meanv,axis=1),axis=1)
		self.bonds = meanv
		return self



	'''
	# Rotate the whole geometry counterclockwise about arbitrary axis
	'''
	def rotate(self,axis,angle):
		rotmat = rotation_matrix(axis,angle)
		self.bonds  = self.bonds @ rotmat
		return self

	'''
	# Taking the index of a bond in the geometry and a target vector, aligns the whole geometry with the target vector
	# Computes the angle between both vectors and rotates the geometry about their cross product
	'''
	def orient(self,refbond_id,target_v):
		refbond = self.bonds[refbond_id]
		ang = -do_ang(refbond,target_v)
		axis = np.cross(refbond,target_v)
		self.rotate(axis,ang)
		return self

	'''
	# The azimuthal Rho is defined as the angle between the u-B0 plane and the v-u plane
	# We can reorient the geometry by rotating around u so that Rho == 0
	# Impossible if u is aligned with B0 (in which case the DC is the same for all values of Rho)
	'''
	def reorientAz(self):
		u = self.bonds[self.tiltID]
		v = self.bonds[self.azID]
		unorm = u/np.linalg.norm(u)
		B0norm = self.B0/np.linalg.norm(self.B0)
		if np.allclose(unorm,B0norm,atol=1e-4):
			warnings.warn("Could not orient Rho as the tilt vector points in the same direction as B0")
			return self
		B0_plane = np.cross(u,self.B0)
		# Invert the definition of B0 plane
		B0_plane = np.cross(self.B0,u)

		v_plane = np.cross(v,u)
		ang = -do_ang(B0_plane,v_plane)
		# This works if v_plane points in the direction of B0, otherwise the rotation must be clockwise
		# MUST be inverted if direction of B0 plane inverted
		if np.dot(self.B0,v_plane) > 0:
			ang *= -1.
		self.rotate(u,ang)
		return self


	'''
	# Take in a list of indices of the bonds that determine the DC, and calculate it
	'''
	def calc_DC(self,fromArray=None,invertSecondOrderTau=True):
		dc = 0.0
		if fromArray is None:
			for i in self.DCbond_IDs:
				bond = self.bonds[i]
				cosang = do_ang(bond,self.B0,returnCos=True)
				dc += self.kDC*0.5*(3*(cosang)**2 - 1.0)
			dc /= len(self.DCbond_IDs)
		# Special case - using secondOrderTau
		elif fromArray == 'tau2':
			for tauVID in self.tilt2ID:
				tauV = self.bonds[tauVID]
				if invertSecondOrderTau:
					tauV = -tauV
				taucos = do_ang(tauV,self.B0,returnCos=True)
				tauidx = find_nearest(self.costau,taucos)
				dc += self.TAU2proj[tauidx]
			dc /= len(self.tilt2ID)
		else:
			taucos = self.calc_cosTau()
			rhocos = self.calc_cosRho()
			tauANG = self.calc_TauANG()
			tauidx = find_nearest(self.costau,taucos)
			rhoidx = find_nearest(self.cosrho,rhocos)
			tauANGidx = find_nearest(self.tauANG,tauANG)
			cases = {'tau-rho':[self.DC_surface,(tauidx,rhoidx)],
					'tau':[self.TAUproj,(tauidx)],
					'rho':[self.RHOproj,(rhoidx)],
					'tauANG-rho':[self.DC_surface,(tauANGidx,rhoidx)]}
			case = cases[fromArray]
			DCarray = case[0]
			DCidx = case[1]
			dc += DCarray[DCidx]
		return dc

	def explore_tiltaz_surface(self,fromArray=None,useCosTau=False,BC_ders=False):
		Ntilt = self.gridsize[0]
		for i in range(Ntilt):
			progress = int(i*100/Ntilt)
			if not progress % 5:
				print("Exploring the Tau-Rho surface...",int(progress),"%", end='\r')
			self.orient(self.tiltID,self.tiltRotv)
			self.reorientAz()
			self.orient(self.tiltID,self.B0)
			if useCosTau:
				dtau = self.tau[i]
			else:
				dtau = self.tauANG[i]
			self.rotate(self.tiltRotv,dtau)
			for j in range(self.gridsize[1]):
				# This assumes that Rho does not affect the DC if u == B0
				if i and i < self.gridsize[0]-1:
					self.reorientAz()
				drho = self.atan2rho[j]
				self.rotate(self.bonds[self.tiltID],drho)
				dc = self.calc_DC(fromArray=fromArray)
				self.DC_surface[i,j] = dc
				if BC_ders:
					ders = self.calc_XYZderivatives()
					for k, d in enumerate(ders):
						self.BC_ders[k,i,j] = d

		self.TAUproj = np.mean(self.DC_surface,axis=1)
		self.RHOproj = np.mean(self.DC_surface,axis=0)
		return self
    
    
	def calc_cosTau(self):
		return do_ang(self.bonds[self.tiltID],self.B0,returnCos=True)

	def calc_TauANG(self):
		return do_ang(self.bonds[self.tiltID],self.B0,returnCos=False)

	def calc_cosRho(self):
		u = np.cross(self.bonds[self.azID],self.bonds[self.tiltID])
		v = np.cross(self.bonds[self.tiltID],self.B0)
		return do_ang(u,v,returnCos=True)

	def calc_atan2Rho(self):
		u = self.bonds[self.tiltID]
		v = self.bonds[self.azID]
		vxu = np.cross(v,u)
		b0xu = np.cross(self.B0,u)
		dotrho = np.dot(vxu,b0xu)
		unorm = u / np.linalg.norm(u)
		detrho = np.dot(unorm,np.cross(vxu,b0xu))
		atan2rho = np.arctan2(detrho,dotrho)
		return atan2rho

	def calc_sineRho(self):
		u = np.cross(self.bonds[self.azID],self.bonds[self.tiltID])
		v = np.cross(self.bonds[self.tiltID],self.B0)
		return do_sine(u,v)

	def eval_grad(self,remove_edges=True,cutoff=0.99,halfRho=False):
		grad = np.gradient(self.DC_surface)
		self.DCGrad_TAU = grad[0]
		self.DCGrad_RHO = grad[1]
		self.DCGrad_TAUproj = np.mean(self.DCGrad_TAU,axis=1)
		self.DCGrad_RHOproj = np.mean(self.DCGrad_RHO,axis=0)
		# Half the derivatives of Rho since the range is twice that of Tau (normalizes them in that case)
		if halfRho:
			self.DCGrad_RHO /= 2.
		if remove_edges:
			cutoff = (cutoff*np.array(self.gridsize)).astype(int)
			for i in range(self.gridsize[0]):
				tauAng = self.tauANG[i]
				for j in range(self.gridsize[1]):
					atan2rho = self.atan2rho[j]
					cutTau = i >= cutoff[0] or i < self.gridsize[0] - cutoff[0]
					cutRho = j >= cutoff[1] or j < self.gridsize[1] - cutoff[1]
					if cutTau or cutRho:
						self.DCGrad_TAU[i,j] = 0.0
						self.DCGrad_RHO[i,j] = 0.0
						if cutTau:
							self.DCGrad_TAUproj[i] = 0.0
						if cutRho:
							self.DCGrad_RHOproj[j] = 0.0
		return self

	def write_array(self,array,suffix="_DC"):
		o = open(self.name+suffix,'w')
		for i in range(self.gridsize[0]):
			o.write("{")
			for j in range(self.gridsize[1]):
				if j:
					o.write(',')
				o.write('%.10E' % Decimal(str(array[i,j])))
			o.write("}")
			if i < self.gridsize[0]-1:
				o.write(",")
			o.write("\n")
		np.savetxt("{}.txt".format(self.name+suffix),array)


