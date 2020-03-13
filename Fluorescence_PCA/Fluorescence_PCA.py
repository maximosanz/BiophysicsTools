import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

class Fluorescence_Data:
    def __init__(self,fn,sep=','):
        self.Spectra, self.Compounds, self.Samples = _parse_file(fn,sep=',')
        self.Spectra_Norm = _normalise_spectra(self.Spectra)
        self.Sorted_Samples = _sorted_unique(self.Samples)
        self.Colours = [ "C{}".format(i) if s.lower() != "buffer" else "gray" for i, s in enumerate(self.Sorted_Samples) ]
        return
    
    def do_PCA(self,
               use_compounds=None,
               use_samples=None,
               exclude_buffer=True,
               N_Components=2,
               normalised=True,
               remove_columns=[]):
        compound_idx = None
        self.pca_compounds = self.Compounds
        self.pca_samples = self.Samples
        if use_compounds is not None:
            compound_idx = [ self.Compounds.index(c) for c in use_compounds ]
            self.pca_compounds = use_compounds
        sp = self.Spectra
        if normalised:
            sp = self.Spectra_Norm
        X = _concat_spectra(sp,compound_idx)
        if remove_columns:
            mask = np.zeros(X.shape[0],dtype=bool)
            mask[remove_columns] = True
            X = X[~mask]
            self.pca_samples = list(np.array(self.pca_samples)[~mask])
        if use_samples is not None:
            mask = np.isin(np.array(self.pca_samples),use_samples)
            self.pca_samples = list(np.array(self.pca_samples)[mask])
            X = X[mask]
        elif exclude_buffer:
            isBuffer = np.array([ True if s.lower() == "buffer" else False for s in self.pca_samples ])
            X = X[~isBuffer]
            self.pca_samples = list(np.array(self.pca_samples)[~isBuffer])
        self.pca = PCA(n_components=N_Components)
        self.pca_proj = self.pca.fit_transform(X)
        return

    def plot_PCA(self,savefig=""):
        fig, axes = plt.subplots(ncols=1)
        Samples = np.array(self.pca_samples)
        for s in _sorted_unique(self.pca_samples):
            IDX = np.where(Samples==s)
            col = self.Colours[self.Sorted_Samples.index(s)]
            plt.scatter(self.pca_proj[IDX,0],self.pca_proj[IDX,1],label=s,color=col,s=70)
        plt.subplots_adjust(right=0.7)
        plt.xlabel("Principal Component 1 ({}%)".format(int(self.pca.explained_variance_ratio_[0]*100)),fontsize=16)
        plt.ylabel("Principal Component 2 ({}%)".format(int(self.pca.explained_variance_ratio_[1]*100)),fontsize=16)
        fig.legend(loc=7,fontsize=14)
        if savefig:
            plt.savefig(savefig)
        plt.show()
        return
        
    def plot_Compound_contribution(self,savefig=''):
        if savefig:
            savefigfn = ".".join(savefig.split('.')[:-1])
            savefigsuf = savefig.split('.')[-1]
        components = self.pca.components_
        N_Components = components.shape[0]
        N_Compounds = len(self.pca_compounds)
        contrib = np.zeros((N_Components,N_Compounds))
        ct = 0
        X_ticks = []
        fig, ax = plt.subplots(N_Components,1)
        for i in range(N_Compounds):
            cmp = self.pca_compounds[i]
            cmp_idx = self.Compounds.index(cmp)
            N_Points = self.Spectra[cmp_idx].shape[1]
            contrib[:,i] = np.absolute(components[:,ct:ct+N_Points]).sum(1)
            for j in range(N_Components):
                if not ct:
                    ax[j].axvline(x=0,color='k')
                ax[j].axvline(x=ct+N_Points,color='k')
            X_ticks.append(ct+N_Points/2.)
            ct += N_Points
        for j in range(N_Components):
            ax[j].bar(np.arange(components.shape[1]),components[j])
            ax[j].set_xticks(X_ticks)
            ax[j].set_xticklabels(self.pca_compounds,fontsize=12)
            ax[j].set_ylabel("PC {}".format(j+1),fontsize=12)
        ax[-1].set_xlabel("Compound",fontsize=12)
        if savefig:
            plt.savefig(savefigfn+'_0.'+savefigsuf)
        plt.show()
        fig, ax = plt.subplots(N_Components,1)
        for j in range(N_Components):
            ax[j].bar(self.pca_compounds,contrib[j])
            ax[j].set_ylabel("PC {}".format(j+1),fontsize=12)
        ax[-1].set_xlabel("Compound",fontsize=12)
        if savefig:
            plt.savefig(savefigfn+'_1.'+savefigsuf)
        plt.show()
        return
        
    def plot_spectra(self,plot_compound,plot_sample=[],plot_columns=[],normalised=True,savefig=''):
        ax = plt.gca()
        color=next(ax._get_lines.prop_cycler)['color']
        sp = self.Spectra
        if normalised:
            sp = self.Spectra_Norm
        Samples = np.array(self.Samples)
        compound_IDX = self.Compounds.index(plot_compound)
        W = sp[compound_IDX][0]
        X = sp[compound_IDX][1:]
        if not plot_sample:
            plot_sample = self.Sorted_Samples
        for s in plot_sample:
            col = self.Colours[self.Sorted_Samples.index(s)]
            IDX = np.where(Samples==s)
            N_lines = X[IDX].shape[0]
            for i in range(N_lines):
                label = None
                if not i:
                    label = s
                plt.plot(W,
                         X[IDX][i],
                         label=label,
                         color=col)
        plt.xlabel("Wavelength (nm)",fontsize=14)
        plt.ylabel("Emission (au)",fontsize=14)
        plt.legend(loc='best',fontsize=12)
        if savefig:
            plt.savefig(savefig)
        plt.show()
        
    def which_column(self,sample):
        for i in range(len(self.Samples)):
            if self.Samples[i] == sample:
                print("{} matches column number {}".format(sample,i))
        
        
def _sorted_unique(l):
    return sorted(set(l), key=l.index)

def _parse_file(fn,sep=','):
    ls = open(fn).readlines()
    header = ls[0]
    Samples = [ s.strip() for s in header.split(sep)[1:] ]
    Compounds = []
    Spectra = []
    current_spectrum = None
    for l in ls[1:]:
        c = l.split(sep)
        if c[0].upper() == "COMPOUND":
            if current_spectrum is not None:
                Spectra.append(np.array(current_spectrum,dtype=np.float64).T)
            current_spectrum = []
            Compounds.append(c[1].strip())
            continue
        current_spectrum.append(c)
    if current_spectrum is not None:
        Spectra.append(np.array(current_spectrum,dtype=np.float64).T)
    return Spectra, Compounds, Samples
        
def _normalise_spectra(Spectra):
    Spectra_Norm = []
    for sp in Spectra:
        X = np.copy(sp)
        X = X.T
        X[:,1:] -= X[:,1:].min(0)
        X[:,1:] /= X[:,1:].max(0)
        X = X.T
        Spectra_Norm.append(X)
    return Spectra_Norm
        
def _concat_spectra(Spectra,idx=None):
    if idx is None:
        idx = list(range(len(Spectra)))
    concat_l = [ Spectra[i][1:] for i in idx ]
    return np.concatenate(concat_l,axis=1)
        

