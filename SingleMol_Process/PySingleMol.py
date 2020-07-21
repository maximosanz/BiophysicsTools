import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.signal import convolve
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import matplotlib.pyplot as plt
import pywt


class PhotonData:

	def __init__(	self,
					fn,
					acceptor_colID=0,
					donor_colID=1,
					N_channels=2,
					donor_color='green',
					acceptor_color='magenta'):
		self.X_raw = np.array([l.split()[:N_channels] for l in open(fn)],dtype=float)
		self.time = np.arange(self.X_raw.shape[0]) / 1000. # In seconds
		self.acceptor_colID = acceptor_colID
		self.donor_colID = donor_colID
		self.baseline = None
		self.N_channels = N_channels
		self.X = np.copy(self.X_raw)
		self.donor_color=donor_color
		self.acceptor_color=acceptor_color
		self.histogram = None
		return

	def old_processing(	self,
						Bleeding_correction=0.05):
		# The following reproduces what the old matlab script does:
		# 1. Background correction by taking the mean:
		X_old = self.X_raw - np.floor(np.mean(self.X_raw,axis=0))
		# 2. Subtracting 5% of the donnor channel from acceptor to account for bleeding
		X_old[:,self.acceptor_colID] -= X_old[:,self.donor_colID] * Bleeding_correction
		self.X = X_old
		return self

	def compute_baseline(	self,
							MinPooling_windowSize=1500,
							convolution_filterSize=14999):
		# Create a baseline by following the rolling minimum on a determined window size
		rolling_min = -maximum_filter(-self.X_raw,size=(MinPooling_windowSize,1))
		# Convolution operation to smooth the baseline:
		# The convolution_filterSize has to be uneven to keep the symmetry
		baseline = convolve(	rolling_min[:,],
											np.ones((convolution_filterSize,1))/convolution_filterSize,
											mode='same')
		# Correct the edges, which are distorded in the convolution 
		edgeSize = int(convolution_filterSize/2)
		baseline[:edgeSize,:] = baseline[edgeSize,:]
		baseline[-edgeSize:,:] = baseline[-edgeSize,:]
		self.baseline = baseline
		return self

	def baseline_correct(	self,
							mintozero=False):
		if self.baseline is None:
			self.compute_baseline()
		self.X = self.X_raw - self.baseline + np.mean(self.baseline,axis=0)
		if mintozero:
			self.X -= np.min(self.X,axis=0)
		return self

	def calc_histogram(	self,
						bins=200,
						Hist_lims=None,
						density=False):
		self.histogram = []
		for i in range(self.X.shape[1]):
			X = self.X[:,i]
			if Hist_lims is None:
				HRange = (X.min(),X.max())
			else:
				HRange = Hist_lims
			self.histogram.append(np.histogram(X,bins,range=HRange,density=density))
		return self

	def find_peaks( self,
					method='wavelet',
					cutoff=None,
					sigmacutoff=2,
					Nlevels=3,
					wavelet = 'rbio3.5',
					sensitivity=0.04):
		Methods = ['wavelet','cutoff']
		if method not in Methods:
			print("ERROR: No such method '{}' to find peaks. Available methods are:".format(method),Methods)
		elif method == 'cutoff':
			if cutoff is not None:
				eps = cutoff
			else:
				eps = np.mean(self.X,axis=0)+np.std(self.X,axis=0)*sigmacutoff
			self.ispeak = self.X > eps
		elif method == 'wavelet':
			self.ispeak = []
			# Perform a wavelet decomposition on each channel, with the specified sensitivity on the first component
			for col in range(self.X.shape[1]):
				X = self.X[:,col] - np.nanmin(self.X[:,col])
				coeff = pywt.wavedec(X, wavelet,level=Nlevels)
				for i in range(Nlevels+1):
					eps = 1.0
					if not i:
						eps = sensitivity
					eps *= np.nanmax(coeff[i])
					coeff[i] = pywt.threshold(coeff[i], value=eps, mode="soft" )
				# Mark non-zero values as peaks (> 0.1)
				reconstructed_signal = np.absolute(pywt.waverec(coeff, wavelet))
				self.ispeak.append(reconstructed_signal > 0.1)
			self.ispeak = np.array(self.ispeak).T
		return self

	def count_peaks(self,threshold=10):
		# This function separates the individual peaks into lists
		# Groups together peaks that are closer than threshold (ms) to each other (to correct artifacts from the wavelet transform)
		self.peaks = []
		self.peak_map = np.empty(self.X.shape)
		self.peak_map[:] = np.nan
		for dim in range(self.X.shape[1]):
			peakIdx = np.where(self.ispeak[:,dim])[0]
			self.peaks.append([])
			Npeaks = len(peakIdx)
			pk_distMat = np.absolute(peakIdx - np.expand_dims(peakIdx,-1))
			flat_distMat = np.triu(pk_distMat).flatten()
			condensed_distMat = flat_distMat[np.nonzero(flat_distMat)]
			Z = linkage(condensed_distMat)
			labels = fcluster(Z,threshold,'distance')
			for i in range(np.max(labels)):
				myPk = peakIdx[labels==i+1]
				pkMin = np.min(myPk)
				pkMax = np.max(myPk)
				integral = self.X[pkMin:pkMax,dim].sum()
				self.peaks[dim].append((pkMin,pkMax,integral))
				self.ispeak[pkMin:pkMax,dim] = True
				self.peak_map[pkMin:pkMax,dim] = i
		return self


	def calc_Q(self):
		N_Colocalised = np.all(self.ispeak,axis=1).sum()
		N_Peaks = self.ispeak.sum(0)
		self.Q = (N_Colocalised - (N_Peaks/self.ispeak.shape[0])*N_Peaks)/N_Peaks
		return self

	def plot_X(self,show=True,savefig=None):
		plt.clf()
		plt.plot(self.time,self.X[:,self.donor_colID],color=self.donor_color)
		plt.plot(self.time,-self.X[:,self.acceptor_colID],color=self.acceptor_color)
		high = np.max(self.X)*1.1
		plt.ylim(-high,high)
		plt.ylabel("Number of photons")
		plt.xlabel("Time (s)")
		if savefig is not None:
			plt.tight_layout()
			plt.savefig(savefig)
		if show:
			plt.show()
		return self
		
	def plot_histogram(	self,
						donor=True,
						acceptor=True,
						savefig=None,
						show=True):
		plt.clf()
		if self.histogram is None:
			self.calc_histogram()
		alpha=1.0
		if donor and acceptor:
			alpha = 0.5
		if donor:
			H = self.histogram[self.donor_colID]
			barw = H[1][1]-H[1][0]
			plt.bar((H[1][1:]+H[1][:-1])/2,H[0],width=barw,color=self.donor_color,alpha=alpha,linewidth=0)
		if acceptor:
			H = self.histogram[self.acceptor_colID]
			barw = H[1][1]-H[1][0]
			plt.bar((H[1][1:]+H[1][:-1])/2,H[0],width=barw,color=self.acceptor_color,alpha=alpha,linewidth=0)
		plt.yscale('log')
		plt.xlabel("Number of photons")
		plt.ylabel("Frequency of occurrence")
		if savefig is not None:
			plt.tight_layout()
			plt.savefig(savefig)
		if show:
			plt.show()
		return self

	def plot_peaks( self,show=True,savefig=None):
		plt.clf()
		X_noise = np.copy(self.X)
		X_peaks = np.copy(self.X)
		X_noise[self.ispeak] = np.nan
		X_peaks[~self.ispeak] = np.nan
		plt.plot(self.time,X_noise[:,self.donor_colID],color='gray',alpha=0.2)
		plt.plot(self.time,-X_noise[:,self.acceptor_colID],color='gray',alpha=0.2)
		plt.plot(self.time,X_peaks[:,self.donor_colID],color=self.donor_color)
		plt.plot(self.time,-X_peaks[:,self.acceptor_colID],color=self.acceptor_color)
		high = np.max(self.X)*1.1
		plt.ylim(-high,high)
		plt.ylabel("Number of photons")
		plt.xlabel("Time (s)")
		if savefig is not None:
			plt.tight_layout()
			plt.savefig(savefig)
		if show:
			plt.show()
		return self

	def plot_all_peaks( self):
		X_noise = np.copy(self.X)
		X_peaks = np.copy(self.X)
		X_noise[self.ispeak] = np.nan
		X_peaks[~self.ispeak] = np.nan
		mean_noise = np.nanmean(X_noise,axis=0)
		for dim in range(self.X.shape[1]):
			for pk in self.peaks[dim]:
				print("Peak between positions {} and {} (integral of {})".format(*pk))
				pkwidth = pk[1]-pk[0]
				xlim = (pk[0]-pkwidth*2,pk[1]+pkwidth*2)
				fct = 1.0
				col = self.donor_color
				if not dim:
					fct = -1.0
					col = self.acceptor_color
				plt.fill_between(np.arange(pk[0],pk[1]),np.mean(mean_noise[dim])*fct,self.X[pk[0]:pk[1],dim]*fct,color=col)
				plt.plot(X_noise[:,self.donor_colID],color='gray',alpha=0.2)
				plt.plot(-X_noise[:,self.acceptor_colID],color='gray',alpha=0.2)
				plt.plot(X_peaks[:,self.donor_colID],color=self.donor_color)
				plt.plot(-X_peaks[:,self.acceptor_colID],color=self.acceptor_color)
				plt.xlim(*xlim)
				high = np.max(self.X)*1.1
				plt.ylim(-high,high)
				plt.ylabel("Number of photons")
				plt.xlabel("Time (ms)")
				plt.show()
		return self

	def write_data(self,ofn,include_peaks=False):
		o = open(ofn,'w')
		for i in range(self.X.shape[0]):
			o.write(("{} "*self.X.shape[1]).format(*np.rint(self.X[i]).astype(int)))
			if include_peaks:
				peakstr = np.array(["B"]* self.X.shape[1])
				peakstr[self.ispeak[i]] = "P"
				o.write(("{} "*self.X.shape[1]).format(*peakstr))

			o.write('\n')
		o.close()
		return self

	#def show_statistics(self):
	#	self.calc_Q()
	#	print("Donor Q value = {:4.2f}\nAcceptor Q value = {:4.2f}".format(Data.Q[Data.donor_colID],Data.Q[Data.acceptor_colID]))
	#	return self








