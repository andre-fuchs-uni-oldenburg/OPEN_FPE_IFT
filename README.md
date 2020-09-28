# OPEN_FPE_IFT
This is a user-friendly open-source Matlab package developed by the research group Turbulence, Wind energy and Stochastics (TWiSt) at the Carl von Ossietzky University of Oldenburg (https://uol.de/en/physics/twist). This package enables to perform a standard analysis of given turbulent data and extracts the stochastic equations describing the scale-dependent cascade process in turbulent flows through Fokker-Planck equations. As the analysis of the scale-dependent cascade process through a hierarchy of spatial and temporal scales in turbulent flows is an integral part of turbulence theory, this interdisciplinary treatment of the turbulent cascade process has the potential for a new way to link the statistical description of turbulence (via common two-point increment statistics), non-equilibrium stochastic thermodynamics and local turbulent flow structures. The presented package can be used also for the analysis of other data with turbulent like complexity.

## Peer-reviewed paper on OPEN_FPE_IFT:
Journal of Open Research Software (under review)

## Standalone application (no MATLAB license needed) 
Standalone applications (64-bit) for Windows®, Linux®, and macOS are also created to run the Matlab code on target machines that do not have MATLAB installed.
When running the standalone applications, MATLAB Runtime 2020a is installed automatically. If a problem occurs during the installation of MATLAB Runtime 2020a, see: https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html

## Getting Started (MATLAB license needed)
The software was implemented in MATLAB 2020a, but it has also been tested on MATLAB 2017b. Before using this script the following toolboxes should be included in your MATLAB license.
1. Curve Fitting Toolbox
2. Optimization Toolbox
3. Parallel Computing Toolbox
4. Signal Processing Toolbox
5. Statistics and Machine Learning Toolbox  

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
This program is distributed under the terms of the Creative Commons Attribution 4.0 International License (CC-BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited. See http://creativecommons.org/licenses/by/4.0/.

### In case of difficulties
Please post an issue on the GitHub repository.
