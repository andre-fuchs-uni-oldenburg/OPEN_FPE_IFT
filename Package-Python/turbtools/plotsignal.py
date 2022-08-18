#!/usr/bin/env python3

"""
Python Package for Statistical Analysis of Turbulence Data
url="https://github.com/aakash30jan/turbtools/"
version="22.6a1"
license="GNU General Public License, version 3"
"""

import numpy as np
from scipy.stats import describe


import matplotlib.pyplot as plt
#config



def plot_signal(data, outpath):
    print("Raw signal description: ", describe(data))
    plt.clf()
    plt.plot(np.arange(data.size), data)
    plt.savefig(outpath+'plot_time_vs_data.png')
    plt.show()
    print('plotsaved to '+outpath+'plot_time_vs_data.png')
    return
