# turbtools
Python Package for Statistical Analysis of Turbulence Data

## Installation
### PyPI [![Downloads](https://pepy.tech/badge/turbtools)](https://pepy.tech/project/turbtools)
`turbtools` is available as a python package from [https://pypi.org/project/turbtools/](https://pypi.org/project/turbtools/).
Download and install it as a system or environment package with pip. It can be then used in a CLI mode or as a python module 
```bash
$ pip install turbtools
```

### Source
Alternatively, the latest `turbtools` package source tarball can be downloaded from [here](https://github.com/aakash30jan/turbtools/archive/refs/heads/main.zip) (18.9 kB). 

## Usage

```console
Usage: turbtools [-h] [-o] [-f datafile] FUNCTION_NAME

Options
  -h, --help            Show this help message and exit
  -o, --outpath         Path to save output stats and plots
  -f, --datafile        Path of datafile in .mat or .nc
  FUNCTION_NAME         Name of the the function
  
```

### Examples
Initialize the datafile to analyze
```console
$ wget https://github.com/andre-fuchs-uni-oldenburg/OPEN_FPE_IFT/raw/master/Renner_8000_Hz.mat
$ turbtools -f ./Renner_8000_Hz.mat
```

Describe stats and plot the raw signal
```console
$ turbtools plot_signal
```

Plot stats stationarity
```console
$ turbtools plot_stationarity
```

Plot probability density function
```console
$ turbtools plot_pdf
```

Plot turbulent energy spectrum
```console
$ turbtools plot_spectrum
```

 
## Issues:
Problems? Please raise an issue at [https://github.com/aakash30jan/turbtools/issues](https://github.com/aakash30jan/turbtools/issues).

[![Issues](https://img.shields.io/github/issues/aakash30jan/turbtools)](#turbtools)  [![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](#turbtools)



## License
This work is licensed under a GNU General Public License Version 3 . [![Open Source Love svg3](https://badges.frapsoft.com/os/v3/open-source.svg?v=103)](#turbtools)



