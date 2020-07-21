#/usr/bin/env python

import PySingleMol
import argparse

DEFAULT_SENSITIVITY = '0.04'

argparser = argparse.ArgumentParser(description='Process Photon-counting Single Molecule data')

argparser.add_argument('-datafile',metavar='PhotonCount_Data',help='Input file with the photon count data',required=True)
argparser.add_argument('-baseline_correction',help='Apply a baseline correction to remove the instrument drift)',action='store_true',default=False)
argparser.add_argument('-pick_signal',help='Pick the signal using the wavelet-based peak detection algorithm',action='store_true',default=False)
argparser.add_argument('-wavelet_sensitivity',metavar='S',help='Sensitivity for the wavelet-based peak detection algorithm (default = {})'.format(DEFAULT_SENSITIVITY),default=DEFAULT_SENSITIVITY)
argparser.add_argument('-separate_histograms',help='Save one histogram for each channel',action='store_true',default=False)
argparser.add_argument('-outdata',metavar='Processed_Data',help='Output data file after processing',default=None)
argparser.add_argument('-outpeaks',metavar='Peaks_Data',help='Output data with detected peaks',default=None)
argparser.add_argument('-outplot',metavar='Peaks_Data',help='Output data with detected peaks',default=None)
argparser.add_argument('-outhisto',metavar='Histogram_Figure',help='Histogram output figure filename',default=None)

args = argparser.parse_args()

Data = PySingleMol.PhotonData(args.datafile)

if args.baseline_correction:
	Data.baseline_correct()

if args.pick_signal or args.outpeaks is not None:
	sensitivity = float(DEFAULT_SENSITIVITY)
	Data.find_peaks(method='wavelet',Nlevels=3,sensitivity=sensitivity)

if args.outdata is not None:
	Data.write_data(args.outdata)

if args.outpeaks is not None:
	Data.write_data(args.outpeaks,include_peaks=True)

if args.outplot is not None:
	if args.pick_signal:
		f = Data.plot_peaks
	else:
		f = Data.plot_X
	f(show=False,savefig=args.outplot)
