#!/usr/bin/env python3

"""
An open source package to perform basic and advanced statistical analysis of turbulence data and other complex systems
url="https://github.com/aakash30jan/turbtools/"
version="22.8a1"
license="GNU General Public License, version 3"
"""

from __future__ import print_function
import sys
import argparse
import os
import string
import requests

import numpy as np
import h5py

from scipy.stats import describe
import scipy.io as sio

import matplotlib.pyplot as plt
#config



#cache location
#IDEALLY IN $HOME
turbtools_cache = './.turbtools_cache/'
if not os.path.exists(turbtools_cache):
    try:
        os.mkdir(turbtools_cache)
    except:
        print("> Permission issues. Please ensure you have proper R/W permissions. Stopping.")


def init_dataset(datafile):
    if not os.path.exists(datafile):
        print('> Data file not found at: ', datafile)
        sys.exit()
    else:
        print('> Initializing data file ', datafile )
        ff = open(turbtools_cache+'DATAFILENAME','w')
        ff.write(datafile)
        ff.close()
    return;
    
def load_dataset(datafile):
    print('> Loading dataset from ', datafile)
    try:
        data = sio.loadmat(datafile)
    except:
        arrays = {}
        f = h5py.File(datafile)
        for k, v in f.items():
            arrays[k] = np.array(v)
        data = arrays['data']  
    print('> Loaded data.shape : ', data.shape)
    ##CLEAN AND RESHAPE
    data = data.squeeze()
    return data;

def plot_signal(data, outpath):
    print("> Raw signal description: \n\t", describe(data))
    plt.clf()
    plt.plot(np.arange(data.size), data)
    plt.xlabel(r'Samples in time')
    plt.ylabel(r'Signal intensity')
    plt.savefig(outpath+'plot_time_vs_data.png')
    plt.show()
    print('>> Plot saved to '+outpath+'plot_time_vs_data.png')
    return


from scipy.interpolate import interp1d
def fillmissing(data, method='linear'):
    """
    method = has to be one of 'linear', 'nearest', 'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', or 'next'. 
    'zero', 'slinear', 'quadratic' and 'cubic' refer to a spline interpolation of zeroth, first, second or third order
    'previous' and 'next' simply return the previous or next value of the point 
    'nearest-up' and 'nearest' differ when interpolating half-integers in that 'nearest-up' rounds up and 'nearest' rounds down.

    Example:
    ARR = np.arange(20.0)
    ARR[0] = np.nan; ARR[15] = np.nan
    fillmissing(ARR, method='linear')
    """
    aindexes = np.arange(data.shape[0])
    agood_indexes, = np.where(np.isfinite(data))
    f = interp1d(agood_indexes, data[agood_indexes], bounds_error=False, copy=False, fill_value="extrapolate", kind=method)
    return f(aindexes)

from plot_stationarity import *

def main():
    parser = argparse.ArgumentParser(description='turbtools - Python Package for Statistical Analysis of Turbulence and other Complex Systems Data',epilog='Cheers :)') 
    parser.add_argument("funcname", type=str, default='plot_signal', nargs='?', help="Type name of the analysis of interest. Available: \n\t plot_signal \n\t plot_stationarity  \n\t plot_pdf  \n\t plot_spectrum  \n\t " )
    parser.add_argument("-f", "--datafile", type=str, default='NONE' , required=False, help="Set datafile to be loaded for analysis" )
    parser.add_argument("-o", "--outpath", type=str, default='NONE' , required=False, help="Output path for plots and stats" )

    args = parser.parse_args()
    #Take care of parsed arguments
    funcname = args.funcname
    datafile = args.datafile
    outpath = args.outpath

    #INIT datafile only once
    if datafile != 'NONE':
        init_dataset(datafile)

    #LOAD datafile from cached location
    try:
        datafilename = turbtools_cache+'DATAFILENAME'
        datafile = open(datafilename).readlines()[0]
        data = load_dataset(datafile)
    except:
        print('>> Please initialize with datafile as\n\n\t turbtools -f <DATFILE_WITHPATH>')
        sys.exit()

    #OUTPATH for saving plots and other
    if outpath == 'NONE':
        outpath = './turbtools_results/'
        if not os.path.exists(outpath):
            os.mkdir(outpath)
    print("> Figures and files would be saved to: " , outpath)


    if np.count_nonzero(np.isnan(data)) > 0:
        data = fillmissing(data, method=fillmissingMethod)

    #FUNCNAME    
    if funcname == 'plot_signal':
        plot_signal(data,outpath)
    elif funcname == 'plot_stationarity':
        plot_stationarity(data,outpath,percent=5)
    else:
        pass
        #print('Not implemented: ', funcname)


    if np.count_nonzero(np.isnan(data)) > 0:
        data = fillmissing(data, method=fillmissingMethod)

    # Estimation of the number of bin's
    increment_bin = int(np.ceil(((data.max() - data.min())/np.nanstd(data)*10)))
    if increment_bin%2== 0: increment_bin = increment_bin +1

    fsamp = 12000

    if funcname == 'plot_pdf':
        plot_pdf(data,increment_bin,outpath)
    elif funcname == 'plot_spectrum':
        plotSpectrum(data,outpath,fsamp,increment_bin)
    elif funcname == 'plot_signal' or funcname == 'plot_stationarity':
        pass
    else:
        print('Not implemented: ', funcname)


    """
    m_data  = np.nanmean(data)
    if m_data<0:
        inpStr = raw_input("The mean value is smaller than 0 m/s. Should the mean value be set equal to 1 m/s (see readme)? Enter Y/N ")    
        if inpStr == 'Y':
	        m_data = 1.0
    """


    return ; 

if __name__ == "__main__":
    main()
