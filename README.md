# Spectral Diffusion

# See "ExampleRun.m" in the Example folder for step-by-step code.

# Overview
Code in Matlab 2023a for processing spectral diffusion of multi b-value DWI for multi-compartment diffusion MRI of the kidney. 
It takes a stack of dicoms sorted by b-value trace (nx, ny, bvalue), and approximates the inverse Laplace transform as a spectrum of M exponential components with non negative least squares. 
The regularization factor can be adjusted by hand, or can be set to generalized cross-validation.
The spectral peaks can then be sorted, with one example provided though the provided sorting may not be ideal depending on the application. 
It outputs maps and spectra on a voxel-by-voxel basis as matfiles.

Please reach out to [Mira Liu](https://miramliu.com/) at mirabai.liu@mountsinai.org if you have any questions. 
May 29 2025

# Code Base
This code was adapted by Mira M. Liu from a [simulation](https://github.com/JoaoPeriquito/NNLS_computation_of_renal_DWI) to create spectral maps with an image viewer that shows the spectrum per voxel. It also will sort parameters based on signal fraction and diffusion, and creates multi-compartment flow maps.

The multi-b-value decay curve is fit to an unconstrained sum of exponentials calculated with individual nonnegative least squares fit, i.e. an inverse Laplace transform of each signal decay. It based off of an original multi-exponential fitting model from Thorarin Bjarnason and Ross Mitchell (2010). 
Therefore, if you use any version of this code, please check out, and cite appropriately, ["Estimation of multicomponent flow in the kidney with multi-b-value spectral diffusion"](https://onlinelibrary.wiley.com/doi/10.1002/mrm.30644) by Mira Liu et al., ["Continuous diffusion spectrum computation for diffusion-weighted magnetic resonance imaging of the kidney tubule system"](https://doi.org/10.21037/qims-20-1360) by Joao Periquito & Thomas Gladytz et al. as well as ["AnalyzeNNLS: Magnetic resonance multiexponential decay image analysis"](https://doi.org/10.1016/j.jmr.2010.07.008) by Thorarin Bjarnason & Joseph Mitchel.

This repository is meant to be a code base, not professional software; it is meant to be forked and edited as needed with proper citation.

It is set to assume peak sorting into a maximum of three diffusion compartments, but this maximum can be adjusted on a user basis

The peaks are set to be sorted according to proposed kidney diffusion coefficient boundaries. They can be re-sorted as needed after the Spectral Diffusion is run. This can be done using the output spectral maps which is a volume  (nx, ny, M) of the probability distribution of diffusion coefficients for each (x,y) voxel.

The b-values must be updated according to acquisition, and the input file assumes it is images stacked in the format of (nx, ny, b) of Trace averages with b in ascending order (e.g. b=0, b100, b200... b800). See "StackedDicoms.m" for an example of the format.

The number of basis vectors is set to M=300, but the number and range can be adjusted. 

The regularization factor can be set to either be optimized for every voxel with generalized cross-validation, or to be set to a specific $\lambda$ value if SNR is relatively consistent.

# Alternates
If python for advanced DWI/IVIM is of interest, please check out https://github.com/darksim33/Pyneapple by J. Jasse and T. Gladytz.



# MIT License
Copyright 2025, Mira M. Liu
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
