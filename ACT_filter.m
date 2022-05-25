function denoised_img = ACT_filter(noisy_img, noise_var, threshold_setting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ACT_filter (Version 1.00)  May 2022 (*)
%     Denoising filter based on Adaptive Curvelet Thresholding described in 
%     the paper [#]:
%
% [#] N. Eslahi and A. Aghagolzadeh, "Compressive Sensing Image Restoration
%     Using Adaptive Curvelet Thresholding and Nonlocal sparse Regularization",
%     IEEE Trans. Image Process., vol.25, no.7, pp. 3126-3140, Jul. 2016
%     https://doi.org/10.1109/TIP.2016.2562563
%
% (*) Unlike the preliminary version of ACT_filter, dubbed Curve_denoising 
%     and released in August 2015, this version is able to deal with both 
%     additive stationary white and colored (a.k.a. correlated) noise.  
%     Furthermore, different thresholding schemes, comprehensive syntaxes, 
%     and optimized code for faster execution are included here.
%
%
%   FUNCTION INTERFACE: 
%        denoised_img = ACT_filter(noisy_img, noise_var, threshold_setting) 
%
% ________________________________________________________________________________
%   INPUTS:         |  CLASS:  |  DESCRIPTION:
% --------------------------------------------------------------------------------
%   noisy_img       | (double) |  The noisy image.
% --------------------------------------------------------------------------------
%   noise_var       | (double) |  The noise variance if it is modeled as AWGN, 
%                   |          |  or the noise FFT-PSD (as size as the noisy image).
%                   |          |  If noise_var is not given or set empty, the
%                   |          |  filter will then assume the noise to be AWGN
%                   |          |  and estimate the noise variance internally.
% --------------------------------------------------------------------------------
% threshold_setting |  (char)  |  The selected threshold setting for filtering.
%                   |    OR    |  <options> : <description>
%                   | (string) |     's'  : soft-shrinkage with adaptive threshold
%                   |          |     'h'  : hard-shrinkage with adaptive threshold
%                   |          |  'ksigma': hard-shrinkage with k-sigma threshold
%                   |          |  If threshold_setting is not given or set empty, 
%                   |          |  it will be set to 's'.
%
%
% ________________________________________________________________________________
%   OUTPUT:         |  CLASS:  |  DESCRIPTION:
% --------------------------------------------------------------------------------
%   denoised_img    | (double) |  The denoised image.
%
%
%
%  List of acronyms used in the syntaxes and comments
%    <Acronym>: <Backronym>
%        ACT  : Adaptive Curvelet Thresholding
%        DCuT : Discrete Curvelet Transform
%        AWGN : Additive White Gaussian Noise
%        ML   : Maximum Likelihood
%        MAD  : Median Absolute Deviation
%        PSD  : Power Spectral Density
%   root-PSD  : square root of PSD
%        var  : variance
%        std  : standard deviation
%        coeff: coefficient
%
%
% DEPENDENCY SUBFUNCTIONS:
%     <a href="matlab:help ACT_filter>ACT">ACT</a>
%     <a href="matlab:help ACT_filter>ML_estimator">ML_estimator</a>
%     <a href="matlab:help ACT_filter>cmpt_DCuT_rootPSD">cmpt_DCuT_rootPSD</a>
% 3rd party dependency functions (downloaded from <a href="www.curvelet.org">www.curvelet.org</a>):
%     <a href="matlab:help fdct_wrapping">fdct_wrapping</a>
%     <a href="matlab:help ifdct_wrapping">ifdct_wrapping</a>
%
%
%   AUTHORS: 
%   Nasser Eslahi (nasser.eslahi@tuni.fi, nasser.eslahi@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Example 1: ACT for AWGN
%
%   x = im2double(imread('cameraman.tif'));
%   std_noise  = .2;
%   normalizer = @(n) (n-mean(n(:)))./std(n(:));
%   z = x + std_noise .* normalizer(randn(size(x)));
%   psnr_noisy   = 10*log10((max(x(:))).^2/mean((x(:)-z(:)).^2))
%
%   x_est = ACT_filter(z);
%%%%x_est = ACT_filter(z, std_noise^2);
%   psnr_denoised = 10*log10((max(x(:))).^2/mean((x(:)-x_est(:)).^2))
%
%   figure, subplot(1,3,1), imagesc(x), colormap('gray'),
%   colorbar, title('ground-truth')
%           subplot(1,3,2), imagesc(z), colormap('gray'),
%   colorbar, title(sprintf('noisy image\nPSNR=%0.2fdB',psnr_noisy))
%           subplot(1,3,3), imagesc(x_est), colormap('gray'),
%   colorbar, title(sprintf('denoised image\nPSNR=%0.2fdB',psnr_denoised))
%
% --------------------------------------------------------------------------------
%  Example 2: ACT for stationary correlated noise
%
%   x = im2double(imread('cameraman.tif'));
%   kernel  = .05*rand(7);
%   FFT_PSD = abs(fft2(kernel,size(x,1),size(x,2)).^2)*numel(x);
%   normalizer = @(n) (n-mean(n(:)))./std(n(:));
%   noise = convn(padarray(normalizer(randn(size(x))),(size(kernel)-1)/2,'both','circular'),kernel,'valid');
%   z = x + noise;
%   psnr_noisy   = 10*log10((max(x(:))).^2/mean((x(:)-z(:)).^2))
%
%   x_est = ACT_filter(z, FFT_PSD);
%   psnr_denoised = 10*log10((max(x(:))).^2/mean((x(:)-x_est(:)).^2))
%
%   figure, subplot(1,5,1), imagesc(x), colormap('gray'),
%   colorbar, title('ground-truth')
%           subplot(1,5,2), imagesc(noise), colormap('gray'),
%   colorbar, title('noise')
%           subplot(1,5,3), imagesc(fftshift(FFT_PSD)), colormap('jet'),
%   colorbar, title('global FFT-PSD')
%           subplot(1,5,4), imagesc(z), colormap('gray'),
%   colorbar, title(sprintf('noisy image\nPSNR=%0.2fdB',psnr_noisy))
%           subplot(1,5,5), imagesc(x_est), colormap('gray'),
%   colorbar, title(sprintf('denoised image\nPSNR=%0.2fdB',psnr_denoised))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2022 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only be used for nonprofit 
% noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes is prohibited.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('noise_var','var')||isempty(noise_var)
    est_noise_var = true;  % blind AWGN denoising
else
    est_noise_var = false; % non-blind denoising
end

if ~exist('threshold_setting','var')||isempty(threshold_setting)
    threshold_setting = 's';
end

if ~sum(strcmpi(threshold_setting,{'s','h','ksigma'}))
   error(['Your selected threshold setting ''%s'' does not exist in this code!\n' ...
          'Please check the <a href="matlab:help %s">syntax</a>'], threshold_setting, mfilename);
end

is_real = 1; % consider only real-valued noisy images


%%% forward curvelet transform 
noisy_DCuT_coeffs = fdct_wrapping(noisy_img, is_real);

if est_noise_var
    %%% Estimating noise std by applying the MAD over the high-passed ...
    %%% noisy image to discard the influence of outliers.
    highest_freq_coeff = noisy_DCuT_coeffs{length(noisy_DCuT_coeffs)}{1};
    noise_std = median(abs(highest_freq_coeff(:) - median(highest_freq_coeff(:))))/.6745;  % MAD
    noise_var = noise_std^2;   
end


%%% global FFT-PSD of the noise 
if numel(noise_var) == 1 %% uncorrelated (a.k.a. white) Gaussian noise %%
    FFT_PSD  = noise_var.*ones(size(noisy_img)).*numel(noisy_img);  % flat FFT-PSD of white noise
    noise_type = 'white';
    
else                     %% correlated (a.k.a. colored) Gaussian noise %%
    FFT_PSD  = noise_var; % FFT-PSD of the noise
    if sum(size(FFT_PSD) == size(noisy_img)) ~= 2
        error('The input noise FFT-PSD ''noise_var'' should be of size %dx%d',...
               size(noisy_img,1), size(noisy_img,2))
    end
    if (max(FFT_PSD(:))-min(FFT_PSD(:)))<.015*numel(FFT_PSD)
        noise_type = 'white';
    else
        noise_type = 'colored';
    end
end


%%% computing the noise root-PSD in DCuT domain (or the std of the noise ...
%%% within each subband) given the noise FFT-PSD (or given the noise ...
%%% variance if it is modeled as AWGN).
noise_DCuT_rootPSD = cmpt_DCuT_rootPSD(FFT_PSD);


%%%% Adaptive Curvelet Thresholding 
denoised_DCuT_coeffs = ACT(noisy_DCuT_coeffs, noise_DCuT_rootPSD, threshold_setting, noise_type);


%%% inverse curvelet transform 
denoised_img = ifdct_wrapping(denoised_DCuT_coeffs, is_real);

return









function denoised_DCuT_coeffs = ACT(noisy_DCuT_coeffs, noise_DCuT_rootPSD, ...
                                    threshold_setting, noise_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   *** ACT (Version 2.00)  May 2022
%   Adaptive Curvelet Thresholding described in the paper [*]:
%
%     N. Eslahi and A. Aghagolzadeh, "Compressive Sensing Image Restoration
%     Using Adaptive Curvelet Thresholding and Nonlocal sparse Regularization",
%     IEEE Trans. Image Process., vol.25, no.7, pp. 3126-3140, Jul. 2016
%     https://doi.org/10.1109/TIP.2016.2562563
%
%   FUNCTION INTERFACE: 
%         denoised_DCuT_coeffs = ACT(noisy_DCuT_coeffs, noise_DCuT_rootPSD, ...
%                                    threshold_setting, noise_type)
%
% ________________________________________________________________________________
%   INPUTS:              |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%   noisy_DCuT_coeffs    |  (cell)  | The discrete curvelet coefficients
%                        |          | of the noisy signal
% --------------------------------------------------------------------------------
%   noise_DCuT_rootPSD   |  (cell)  | The noise root-PSD in DCuT domain
%                        |          | (or the std of noise in each subband) 
% --------------------------------------------------------------------------------
%   threshold_setting    |  (char)  | The selected threshold setting for filtering.
%                        |    OR    | <options> : <description>
%                        | (string) |      's'  : ACT via soft-shrinkage
%                        |          |      'h'  : ACT via hard-shrinkage
%                        |          | 'ksigma': hard-shrinkage with k-sigma threshold
% --------------------------------------------------------------------------------
%   noise_type           |  (char)  | The type of noise correlation
%                        |    OR    | 'white'   : uncorrelated Gaussian noise
%                        | (string) | 'colored' : correlated Gaussian noise
%
%
% ________________________________________________________________________________
%   OUTPUT:              |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%   denoised_DCuT_coeffs |  (cell)  | the denoised coefficients
%
%
%   *** ACT (version 1.00) was released in August 2015. In the new version  
%   of ACT, different thresholding schemes, correlated noise removal, and
%   comprehensive syntax are included.
%
%
%   AUTHOR: 
%   Nasser Eslahi (nasser.eslahi@tuni.fi, nasser.eslahi@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% computing the noise DCuT-root-PSD (i.e. the std of noise in each DCuT subband) ...
%%% through scaling the norm2 of the DCuT basis/frame by the noise std.

denoised_DCuT_coeffs = noisy_DCuT_coeffs;

for J = 2:length(noisy_DCuT_coeffs)
    for L = 1:length(noisy_DCuT_coeffs{J})
        
        if sum(strcmpi(threshold_setting,{'s','h'})) %%ACT (Eslahi & Aghagolzadeh, 2016)%%
            
            %%% estimating the sample std of the clean (noise-free) DCuT coefficients ...
            %%% via the maximum likelihood estimator given the noisy DCuT coefficients ...
            %%% and the noise root-PSD in DCuT domain (see Eq.(12) of [*]).
            clean_subband_std = ML_estimator(noisy_DCuT_coeffs{J}{L}, ...
                                             noise_DCuT_rootPSD{J}(L),...
                                             noise_type);


            %%% computing the adaptive threshold for the selected subband and
            %%% then using it within the selected filter
            Threshold = sqrt(2)*((noise_DCuT_rootPSD{J}(L).^2)./clean_subband_std);

            if     strcmpi(threshold_setting,'s')  %% soft-shrinkage %%
                denoised_DCuT_coeffs{J}{L} = sign(noisy_DCuT_coeffs{J}{L}) .* ...
                                          max(abs(noisy_DCuT_coeffs{J}{L}) - Threshold, 0);

            elseif strcmpi(threshold_setting,'h')  %% hard-shrinkage %%
%                 Threshold = max(sqrt(2*log(numel(Threshold))),1)*Threshold./sqrt(2);
            Threshold = (3 + double(J == length(noisy_DCuT_coeffs))) * ...
                                                Threshold./sqrt(2); 
                denoised_DCuT_coeffs{J}{L} = noisy_DCuT_coeffs{J}{L} .* ...
                                        (abs(noisy_DCuT_coeffs{J}{L}) > Threshold);
                                    
            end
            
        
        else  %% k-sigma thresholding (Starck, Candes, Donoho, 2002) %%
            Threshold = (3 + double(J == length(noisy_DCuT_coeffs))) * ...
                                                noise_DCuT_rootPSD{J}(L);                                             
            denoised_DCuT_coeffs{J}{L} = noisy_DCuT_coeffs{J}{L} .* ...
                                    (abs(noisy_DCuT_coeffs{J}{L}) > Threshold);
        end


    end
end

return;









function Est_sample_std_clean_DCuT_coeff = ML_estimator(noisy_DCuT_coeffs,...
                                               noise_root_PSD, noise_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ML_estimator estimates the sample std of the clean (noise-free) coefficients in the DCuT
% domain at a specific subband (scale and rotation angle) via the ML estimator given
% the noisy DCuT coefficients and the noise std.
%
%
% FUNCTION INTERFACE:
%        Est_sample_std_clean_DCuT_coeff = ML_estimator(noisy_DCuT_coeffs,...
%                                              noise_root_PSD, noise_type)
%
% ________________________________________________________________________________
%   INPUTS:             |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%  noisy_DCuT_coeffs    | (double) | Noisy DCuT coefficients at a specific subband
% --------------------------------------------------------------------------------
%  noise_root_PSD       | (double) | Noise DCuT-root-PSD at a specific subband
% --------------------------------------------------------------------------------
%  noise_type           |  (char)  | The type of noise correlation
%                       |    OR    | 'white'   : uncorrelated Gaussian noise
%                       | (string) | 'colored' : correlated Gaussian noise
%
%
% ________________________________________________________________________________
%   OUTPUT:                        |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%  Est_sample_std_clean_DCuT_coeff | (double) | Estimated sample std of clean DCuT 
%                                  |          | coefficients at the selected subband
%
%
% AUTHOR:
%    Nasser Eslahi   (nasser.eslahi@{gmail.com, tuni.fi})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% convolutional kernel for ML averaging (local averaging window)
if strcmpi(noise_type,'white')
    kernel = ones(7)/48; kernel(4,4) = 0;
else
    kernel = ones(31)/960; kernel(16,16) = 0;
end

%%% estimating sample variance of noisy coefficients as the average of ...
%%% squared coefficients over a local moving window
Est_sample_var_noisy_coeff = convn(...
                               padarray(noisy_DCuT_coeffs.^2, (size(kernel)-1)/2, 'both','circular'),...
                                    kernel, 'valid');
                                
%%% estimating the sample variance of clean coefficients by subtracting ...
%%% the noise std from the estimated sample variance of noisy coefficients
Est_sample_var_clean_coeff = Est_sample_var_noisy_coeff - noise_root_PSD.^2;

%%% imposing non-negativity (because it is PSD/variance!)
Est_sample_var_clean_coeff(  Est_sample_var_clean_coeff <0 ) = 0;

%%% taking the square-root for sample std
Est_sample_std_clean_DCuT_coeff = sqrt(Est_sample_var_clean_coeff);
return;









function noise_DCuT_rootPSD = cmpt_DCuT_rootPSD(FFT_PSD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cmpt_DCuT_rootPSD computes the noise root-PSD in DCuT domain (or the std of noise within 
% each subband) given the noise FFT-PSD (or given the noise variance if it is modeled as AWGN)
%
%
% FUNCTION INTERFACE:
%        noise_DCuT_rootPSD = cmpt_DCuT_rootPSD(FFT_PSD)
%
% ________________________________________________________________________________
%  INPUTS:            |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%  FFT_PSD            | (double) | The noise FFT-PSD (as size as the noisy image).
%                     |          | 
%
%
% ________________________________________________________________________________
%  OUTPUTS:           |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%  noise_DCuT_rootPSD |  (cell)  | The noise root-PSD in DCuT domain (or the std  
%                     |          | of the noise within each subband)
%
%
% AUTHOR:
%    Nasser Eslahi   (nasser.eslahi@{gmail.com, tuni.fi})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


is_real = 1;
 
%%% computing the root-PSD in DCuT domain
kernel_noise = fftshift(ifft2(sqrt(FFT_PSD))); % convolutional kernel of the noise
DCuT_frames  = fdct_wrapping(kernel_noise,is_real); % DCuT frames (for all subbands)
func         = @(x) sqrt(mean(x(:).^2)); % mean squared function
noise_DCuT_rootPSD = cellfun(@(c)cellfun(func, c, 'uni', false), DCuT_frames, 'uni', false);
noise_DCuT_rootPSD = cellfun(@cell2mat, noise_DCuT_rootPSD, 'uni', false);

    
return







%%%%% References %%%%%

% E. Candes, L. Demanet, D. Donoho, and X. Ying, "Fast discrete curvelet  
% transforms," Multiscale Model. Simul., vol. 5, no. 3, pp. 861-899, Sep. 2006.

% D. L. Donoho and J. M. Johnstone, "Ideal spatial adaptation by wavelet  
% shrinkage," Biometrika, vol. 81, no. 3, pp. 425-455, 1994.

% M. K. Mihcak, I. Kozintsev, K. Ramchandran, and P. Moulin, "Low-complexity  
% image denoising based on statistical modeling of wavelet coefficients," 
% IEEE Signal Process. Lett., vol. 6, no. 12, pp. 300-303, Dec. 1999.

% D. L. Donoho, "De-noising by soft-thresholding," IEEE Trans. Inf. Theory, 
% vol. 41, no. 3, pp. 613-627, May 1995.

% F. R. Hampel, "The influence curve and its role in robust estimation,"
% J. Am. Stat. Assoc., vol. 69, no. 346, pp. 383–393, 1974.

% J. L. Starck, E. J. Candes and D. L. Donoho, "The curvelet transform for
% image denoising," IEEE Trans. Image Process., vol.11, no.6, pp. 670-684,
% Jun. 2002.

% N. Eslahi and A. Aghagolzadeh, "Compressive sensing image restoration
% using adaptive curvelet thresholding and nonlocal sparse regularization",
% IEEE Trans. Image Process., vol.25, no.7, pp. 3126-3140, Jul. 2016

  