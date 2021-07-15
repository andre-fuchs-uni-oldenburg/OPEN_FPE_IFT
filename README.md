## OPEN_FPE_IFT     
This is a user-friendly open-source Matlab package developed by the research group Turbulence, Wind energy and Stochastics (TWiSt) at the Carl von Ossietzky University of Oldenburg (https://uol.de/en/physics/twist). This package enables to perform a standard analysis of given turbulent data and extracts the stochastic equations describing the scale-dependent cascade process in turbulent flows through Fokker-Planck equations. As the analysis of the scale-dependent cascade process through a hierarchy of spatial and temporal scales in turbulent flows is an integral part of turbulence theory, this interdisciplinary treatment of the turbulent cascade process has the potential for a new way to link the statistical description of turbulence (via common two-point increment statistics), non-equilibrium stochastic thermodynamics and local turbulent flow structures. The presented package can be used also for the analysis of other data with turbulent like complexity.

## README file on OPEN_FPE_IFT:
A detailed description of the individual functions can be found on arXiv: https://arxiv.org/abs/2106.13042

## Standalone application (no MATLAB license needed) 
To enhance the accessibility, standalone applications (64- bit) for Windows, macOS and Linux are also created to run the MATLAB code on target machines that do not have a MATLAB license. When running the standalone applications, MATLAB Runtime 2020a is installed automatically. If a problem occurs during the installation of MATLAB Runtime 2020a, see: https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html

## Getting Started (MATLAB license needed)
The software is implemented in MATLAB 2020a. Before using this script the following toolboxes should be included in your MATLAB license.
1. Curve Fitting Toolbox
2. Optimization Toolbox
3. Parallel Computing Toolbox
4. Signal Processing Toolbox
5. Statistics and Machine Learning Toolbox  
As these MATLAB toolboxes are essential for the application of the package, the compatibility
to Octave can not be provided. 

### Installation  
1. Download the Zip archive and unzip it into an appropriate directory
2. Add the destination folder to you MATLAB search path by: 
```addpath(folderName)```
3. main.m is the primary script in which all the different functions to carry out various analyses are included. All subroutines are called in a logical order.

### Matlab File Exchange
This package can be also directly downloaded from the Matlab File Exchange:

[![View OPEN_FPE_IFT on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://de.mathworks.com/matlabcentral/fileexchange/80551-open_fpe_ift)

## Exemplary dataset 
Example analyses of the dataset: "Renner_8000_Hz.mat" can be found on peer-reviewed paper on OPEN_FPE_IFT.

### License
This open source MATLAB package is available as free software, under the GNU General Public License (GPL) version 3, which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited. See http://creativecommons.org/licenses/by/3.0/.

### In case of difficulties
Please post an issue on the GitHub repository.
