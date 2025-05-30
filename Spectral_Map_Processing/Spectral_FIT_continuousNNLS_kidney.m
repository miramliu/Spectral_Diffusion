%% Fit curve to unconstrained sum of exponentialswith NNLS
% Mira Liu 2025

function [parameter_map, spectralmap] = Spectral_FIT_continuousNNLS_kidney(b_values, ImageStack, lambda)
%------------------------------------------------------%

%------------------------------------------------------%
    [N_Bvalues, nx, ny] = size(ImageStack);

    %% Generate NNLS space of values, not entirely sure about this part, check with TG?
    ADCBasisSteps = 300; %(??)
    ADCBasis = logspace( log10(5), log10(2200), ADCBasisSteps);
    A = exp( -kron(b_values',1./ADCBasis));
    
    %% create empty arrays to fill
    amplitudes = zeros(ADCBasisSteps,1);
    resnorm = zeros(1);
    resid = zeros(length(b_values),1);
    y_recon = zeros(max(b_values),1);
    resultsPeaks = zeros(6,1); %6 was 9 before? unsure why
    
    %create empty parameter map
    parameter_map = zeros(nx,ny,6);
    spectralmap = zeros(nx, ny, 300);

    
    for i=1:nx
        for j=1:ny
            if (ImageStack(1,i,j) > 100 ) 

                % for normal b values
                SignalInput = squeeze(double(ImageStack(1:N_Bvalues,i,j)/ImageStack(1,i,j))); 

               
                %% try to fit them with NNLS
                if strcmp(lambda, 'cv')
                    [TempAmplitudes, TempResnorm, TempResid ] = simpleCVNNLS(A, SignalInput);
                else

                %% fitting with simple NNLS, with an assumed provided constant regularization paramater of lambda = #b-value/SNR 
                    [TempAmplitudes, TempResnorm, TempResid ] = simpleCVNNLS_curveregularized(A, SignalInput, lambda); %this now also still has the ends regularized
                end

                
                amplitudes(:) = TempAmplitudes';
                resnorm(:) = TempResnorm';
                resid(1:length(TempResid)) = TempResid';
                y_recon(1:size(A,1)) = A * TempAmplitudes;
            
                % check resids
                SSResid = sum(resid.^2);
                SStotal = (length(b_values)-1) * var(SignalInput);
                rsq = 1 - SSResid/SStotal; 


                OutputDiffusionSpectrum = amplitudes;

                %% for plotting
                %{
                plot(b_values,SignalInput)
                pause(1)
                %% output renaming, just to stay consistent with the TG&JP code
                
                semilogx((1./ADCBasis)*1000,OutputDiffusionSpectrum)
                hold on;
                xline(.8), xline(5), xline(50);
                pause(1)
                hold off;
                %plot(OutputDiffusionSpectrum);
                %pause(1)
                Chi = resnorm;
                Resid = resid;
            
                %}

                % fit 
                [GeoMeanRegionADC_1,GeoMeanRegionADC_2,GeoMeanRegionADC_3,GeoMeanRegionADC_4,RegionFraction1,RegionFraction2,RegionFraction3,RegionFraction4 ] = NNLS_result_mod_ML_fourpeaks(OutputDiffusionSpectrum, ADCBasis);
                resultsPeaks(1) = RegionFraction1; %(frac_fast - RegionFraction1)./frac_fast.*100;
                resultsPeaks(2) = RegionFraction2; %(frac_med - RegionFraction2)./frac_med.*100;
                resultsPeaks(3) = RegionFraction3; %(frac_slow - )./frac_slow.*100;
                resultsPeaks(4) = RegionFraction4; %(frac_fibro - )./frac_slow.*100;
                resultsPeaks(5) = GeoMeanRegionADC_1; %(diff_fast - GeoMeanRegionADC_1./1000)./diff_fast.*100;
                resultsPeaks(6) = GeoMeanRegionADC_2; %(diff_med - GeoMeanRegionADC_2./1000)./diff_med.*100;
                resultsPeaks(7) = GeoMeanRegionADC_3; %(diff_slow - GeoMeanRegionADC_3./1000)./diff_slow.*100;
                resultsPeaks(8) = GeoMeanRegionADC_4; %(diff_fibro - GeoMeanRegionADC_3./1000)./diff_slow.*100;



                if rsq>0.7
                    if resultsPeaks(1)<1000 %it's set to 10000 if no peaks found, see line 32 of NNLS_result_mod
                        % now  try to sort them... 
                        SortedresultsPeaks = ReSort_fourpeaks(resultsPeaks);
                        parameter_map(i,j,1) = SortedresultsPeaks(1); %vasc frac
                        parameter_map(i,j,2) = SortedresultsPeaks(2); %tubule frac
                        parameter_map(i,j,3) = SortedresultsPeaks(3); %tissue frac
                        % SortedresultsPeaks(4) discard peaks with D < 0.8 (labeled 'fibro')
                        parameter_map(i,j,4) = SortedresultsPeaks(5); %vasc D
                        parameter_map(i,j,5) = SortedresultsPeaks(6); %tubule D
                        parameter_map(i,j,6) = SortedresultsPeaks(7); %tissue D
                        % SortedresultsPeaks(8) discard peaks with D < 0.8 (labeled 'fibro')
    
                        spectralmap(i,j,:) = OutputDiffusionSpectrum;
                    else
                        parameter_map(i,j,:) = zeros(6,1);
                        spectralmap(i,j,:) = zeros(300,1);
                    end
                else
                    parameter_map(i,j,:) = zeros(6,1);
                    spectralmap(i,j,:) = zeros(300,1);
                end

            else
               parameter_map(i,j,:) = zeros(6,1);
               spectralmap(i,j,:) = zeros(300,1);
            end
            
        end
    end
end
