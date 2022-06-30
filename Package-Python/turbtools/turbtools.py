#!/usr/bin/env python3

"""
Python Package for Statistical Analysis of Turbulence Data
url="https://github.com/aakash30jan/turbtools/"
version="22.6a1"
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
    os.mkdir(turbtools_cache)



def init_dataset(datafile):
    if not os.path.exists(datafile):
        print('Data file not found at: ', datafile)
        sys.exit()
    else:
        print('Initializing data file ', datafile )
        ff = open(turbtools_cache+'DATAFILENAME','w')
        ff.write(datafile)
        ff.close()
    return;
    
def load_dataset(datafile):
    print('Loading dataset from ', datafile)
    try:
        data = sio.loadmat(datafile)
    except:
        arrays = {}
        f = h5py.File(datafile)
        for k, v in f.items():
            arrays[k] = np.array(v)
        data = arrays['data']  
    print('Loaded data.shape : ', data.shape)
    ##CLEAN AND RESHAPE
    data = data.squeeze()
    return data;

def plot_signal(data, outpath):
    print("Raw signal description: ", describe(data))
    plt.clf()
    plt.plot(np.arange(data.size), data)
    plt.savefig(outpath+'plot_time_vs_data.png')
    plt.show()
    print('plotsaved to '+outpath+'plot_time_vs_data.png')
    return

def plot_stationarity(data, outpath):
    return 0;

def main():
    parser = argparse.ArgumentParser(description='turbtools - Python Package for Statistical Analysis of Turbulence and other Complex Systems Data',epilog='Cheers :)') 
    parser.add_argument("funcname", type=str, default='plot_signal', nargs='?', help="Type name of the analysis of interest" )
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
        print('Please initialize with datafile as turbtools -f <DATFILE_WITHPATH>')
        sys.exit()

    #OUTPATH for saving plots and other
    if outpath == 'NONE':
        outpath = './turbtools_results/'
        if not os.path.exists(outpath):
            os.mkdir(outpath)

    #FUNCNAME    
    if funcname == 'plot_signal':
        plot_signal(data,outpath)
    elif funcname == 'plot_stationarity':
        plot_stationarity(data,outpath)
    else:
        print('Not implemented: ', funcname)

    return ; 

if __name__ == "__main__":
    main()
