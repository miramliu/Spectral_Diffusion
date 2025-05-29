# Spectral Diffusion

# Overview
Code in Matlab 2023a for processing spectral diffusion of multi b-value DWI for multi-compartment diffusion MRI of the kidney. 
It takes a stack of dicoms sorted by b-value trace (nx, ny, bvalue), and approximates the inverse Laplace transform as a spectrum of M exponential components with non negative least squares. 
The regularization factor can be adjusted by hand, or can be set to generalized cross-validation.
The spectral peaks can then be sorted, with one example provided though the provided sorting may not be ideal depending on the application. 
It outputs maps and spectra on a voxel-by-voxel basis as matfiles.

Please reach out to Mira Liu at mirabai.liu@mountsinai.org if you have any questions.
- May 2025

See "ExampleRun.m" in the Example folder for how to run the code.

# Code Base
This code was adapted by Mira M. Liu from a simulation with [open source code](https://github.com/JoaoPeriquito/NNLS_computation_of_renal_DWI) to create spectral maps with an image viewer that shows the spectrum per voxel. It also will sort parameters based on signal fraction and diffusion, and creates multi-compartment flow maps.

The inverse laplace is calculated with individual nonnegative least squares fit of each signal decay is based off of an original multi-exponential fitting model from Thorarin Bjarnason and Ross Mitchell (2010). 
Therefore, if you use any version of this code, please check out, and cite appropriately, ["Quantification of multi-compartment flow from spectral diffusion of IVIM"](https://doi.org/10.48550/arXiv.2408.06427) by Mira Liu & Jon Dyke et al., ["Continuous diffusion spectrum computation for diffusion-weighted magnetic resonance imaging of the kidney tubule system"](https://doi.org/10.21037/qims-20-1360) by Joao Periquito & Thomas Gladytz et al. as well as ["AnalyzeNNLS: Magnetic resonance multiexponential decay image analysis"](https://doi.org/10.1016/j.jmr.2010.07.008) by Thorarin Bjarnason & Joseph Mitchel.

This repository is meant to be a code base, not professional software; it is meant to be forked and edited as needed with proper citation.

It is set to assume a sort into a maximum of three diffusion compartments, but this maximum can be adjusted on a user basis

The peaks are set to be sorted according to proposed kidney diffusion coefficient boundaries. They can be re-sorted as needed after the Spectral Diffusion is run from the output Spectral maps (SpectralDWI.Spectral_Volume) has an (nx, ny, M) volume of the spectra for each voxel.

The b-values must be updated according to acquisition, and the input file assumes it is images stacked in the format of (nx, ny, b) of Trace averages with b in ascending order (e.g. b=0, b100, b200... b800).

The number of basis vectors is set to 300, but the number and range can be adjusted. 

The regularization factor can be set to either be optimized for every voxel with generalized cross-validation, or to be set to a specific $\lambda$ value if SNR is relatively consistent.

# Alternates
If python for advanced DWI/IVIM is of interest, please check out https://github.com/darksim33/Pyneapple by J. Jasse and T. Gladytz.
