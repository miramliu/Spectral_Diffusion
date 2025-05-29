%% example code to run spectral diffusion of multi-b-value  DWI of a kidney


% load relevant folders, adjust the path to local drive as needed

addpath /Spectral_Diffusion/rNNLS
addpath /Spectral_Diffusion/noise_matrix
addpath /Spectral_Diffusion/Spectral_Map_Processing/
addpath /Spectral_Diffusion/rNNLS/nwayToolbox



% run voxel-by-voxel fit on example

%give total path to folder of interest (folder that has the stacked trace dicoms)
dicompath = '/Spectral_Diffusion/Example/';
StackName = 'StackedDicoms.mat'; % given a stacked set of trace dicoms (nx, ny, b) 
lambda = 0.1; %set to 0.1, can also be set to 'cv' for generalized cross validation, or to alternate values.
KidneyMask = 'KidneyMask.mat';
RunKidney_Spectral_Diffusion(dicompath,StackName,lambda,KidneyMask);

% upon completion should see "Name_SpectralDWI.m" file. 
% structure should contain "Parameter_Volume" which has the three compartment fractions and three compartment diffusion coefficients
% "Spectral_Volume" should contain the spectra per voxel. This can be used to view the spectra per voxel in a GUI that I have in progress. 


% once completed, load the created file (which is saved in the original folder with name_SpectralDWI) 
% and see the fD maps as follows using another GUI of mine: 
load('StackedDicoms_SpectralDWI.mat')
Parameter_Volume = SpectralDWI.Parameter_Volume;
fD_maps(:,:,1)=Parameter_Volume(:,:,1).*Parameter_Volume(:,:,4);
fD_maps(:,:,2)=Parameter_Volume(:,:,2).*Parameter_Volume(:,:,5);
fD_maps(:,:,3)=Parameter_Volume(:,:,3).*Parameter_Volume(:,:,6);
imagestack(fD_maps)


