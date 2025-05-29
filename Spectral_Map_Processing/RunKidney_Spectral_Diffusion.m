% given a path to a mat file, run spectral diffusion analysis and return the parameter map and the map of the spectrum
% Mira Liu 2025

function RunKidney_Spectral_Diffusion(varargin)

    dicomfolderpath = varargin{1};
    StackName = varargin{2};

     % the b-values for the acquisition.
    Bvalues = [0, 10, 30, 50, 80, 120, 200, 400, 800 ];

    loadDicoms = load(fullfile(dicomfolderpath,StackName));
    StackedDicoms = loadDicoms.StackedDicoms;
    
    if nargin > 2
        lambda = varargin{3}; %set lambda to a specific value
        if nargin > 3
            KidneyMaskName = varargin{4}; % give an allograft mask in the same folder to process the kidney only (reduce time)
            KidneyMask = load(fullfile(dicomfolderpath,KidneyMaskName), 'Mask').Mask;
        else
             KidneyMask = ones(size(StackedDicoms(:,:,1))); % no mask
        end
    else
        lambda = 'cv'; %cross validation
        KidneyMask = ones(size(StackedDicoms(:,:,1))); % no mask

    end

    KidneyMaskedDicoms = squeeze(StackedDicoms(:,:,:)).*KidneyMask;
    KidneyMaskedDicoms = permute(KidneyMaskedDicoms, [3,1,2]); %to have it bval, nx, ny

    % run analysis
    disp(['started: '  + string(datetime("now"))])

    [parameter_map, spectral_map] =  Spectral_FIT_continuousNNLS_kidney(Bvalues,KidneyMaskedDicoms,lambda);

    SpectralDWI.Parameter_Volume     = parameter_map;
    SpectralDWI.Spectral_Volume      = spectral_map; %% file size too large... can only save peaks themselves for one slice at a time. Matlab limitation?
    

    % export, can adjust naming preferences here.
    SaveDIR = fullfile(dicomfolderpath, StackName(1:end-4) + "_SpectralDWI.mat"); 
    save (SaveDIR, 'SpectralDWI');
    disp(['saved.... ' SaveDIR])

    disp(['Completed: ' + string(datetime("now"))])
end




