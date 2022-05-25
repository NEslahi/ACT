%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The demo code for filtering of additive stationary Gaussian noise using 
%  the adaptive curvelet thresholding method, proposed in:
%
%     N. Eslahi and A. Aghagolzadeh, "Compressive Sensing Image Restoration
%     Using Adaptive Curvelet Thresholding and Nonlocal sparse Regularization",
%     IEEE Trans. Image Process., vol.25, no.7, pp. 3126-3140, Jul. 2016
%     https://doi.org/10.1109/TIP.2016.2562563
%
%
%   AUTHOR: 
%   Nasser Eslahi (nasser.eslahi@tuni.fi, nasser.eslahi@gmail.com)
%
%                          Release ver. 1.0  (May 2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2022 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only be used for nonprofit 
% noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes is prohibited.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath(genpath('./'));
clear

%%% finding all images in the Test_Images folder
srcFiles   = dir('./Test_Images/*.tif');


%%% loading randomly selected test image
x = im2double(imread(srcFiles(randi([1 numel(srcFiles)])).name)); 
if size(x,3)~=1
   x = rgb2gray(x);
end


%%% generating a stationary Gaussian noise
[noise, noise_FFT_PSD] = generate_Gaussian_noise(size(x));


%%% creating the noisy image
z          = x + noise; 
psnr_noisy = 10*log10((max(x(:))).^2/mean((x(:)-z(:)).^2));


%%% denoising via k-sigma thresholding (Starck, Candes, and Donoho, 2002)
x_est_ksigma = ACT_filter(z, noise_FFT_PSD, 'ksigma');
psnr_ksigma  = 10*log10((max(x(:))).^2/mean((x(:)-x_est_ksigma(:)).^2));


%%% denoising via the adaptive curvelet thresholding (Eslahi and Aghagolzadeh, 2016)
x_est_ACT = ACT_filter(z, noise_FFT_PSD);
psnr_act  = 10*log10((max(x(:))).^2/mean((x(:)-x_est_ACT(:)).^2));


%%% displaying the images
figure, 

subplot(2,3,1), imagesc(noise), colormap('gray'),
colorbar, title('noise'), axis off

subplot(2,3,4), imagesc(fftshift(noise_FFT_PSD)), colormap('jet'),
colorbar, title('noise FFT-PSD'), axis off

subplot(2,3,2), imagesc(x), colormap('gray'),
colorbar, title('ground-truth'), axis off

subplot(2,3,3), imagesc(z), colormap('gray'),
colorbar, axis off, title(sprintf('noisy image\nPSNR=%0.2fdB',psnr_noisy))

subplot(2,3,5), imagesc(x_est_ksigma), colormap('gray'),
colorbar, axis off, title(sprintf('denoised image using k-sigma\nPSNR=%0.2fdB',psnr_ksigma))

subplot(2,3,6), imagesc(x_est_ACT), colormap('gray'),
colorbar, axis off, title(sprintf('denoised image using ACT\nPSNR=%0.2fdB',psnr_act))


%%  Auxiliary function for generating Gaussian noise

function  [noise, noise_FFT_PSD] = generate_Gaussian_noise(SizeX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_Gaussian_noise creates stationary white/colored Gaussian noise
% with respect to a convolutional kernel selected randomly.
%
%
% FUNCTION INTERFACE:
%         [noise, noise_FFT_PSD] = generate_Gaussian_noise(SizeX)
%
% ________________________________________________________________________________
%  INPUT:        |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%  SizeX         | (double) | Size of the ground-truth image.
%
%
% ________________________________________________________________________________
%  OUTPUTS:      |  CLASS:  | DESCRIPTION:
% --------------------------------------------------------------------------------
%  noise         | (double) | generated stationary Guassian noise.
% --------------------------------------------------------------------------------
%  noise_FFT_PSD | (double) | the noise FFT-PSD (as size as SizeX).  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kernel_type  = randi([1 14]); % the randomly selected kernel for noise generation
normalizer   = @(n) (n-mean(n(:)))./std(n(:)); % making the noise zero-mean with unit variance

if kernel_type == 1
    k        = 1:randi([2 11]);
    mu       = mean(k);
    sigma    = rand;
    Gaus     = (1/sqrt(2*pi*sigma^2))*exp(-(k-mu).^2/(2*sigma^2));
    kernel   = Gaus'*Gaus; % 2D separable Gaussian kernel
    kernel   = kernel./norm(kernel(:)); 
  
elseif ismember(kernel_type, 2:3)
    if kernel_type == 2
      kernel = randn(randi([2 11])); % 2D random kernel of Gaussian distribution
    else
      kernel = rand(randi([2 11]));  % 2D random kernel of uniform distribution  
    end
    kernel   = kernel./norm(kernel(:));
       
elseif ismember(kernel_type, 4:5)
    k        = 2*randi([4 12],1)+1;
    kernel   = ceil(k/2) - abs((1:k)-ceil(k/2)); % 1D horizontal kernel
    if kernel_type == 5
      kernel = kernel';    % 1D vertical kernel
    end
    kernel   = kernel / norm(kernel(:)); 
    
elseif ismember(kernel_type, 6:8)
    sz       = randi([16 40],1) * ones(1,2);
    sz2      = -(1 - mod(sz, 2)) * 1 + floor(sz/2);
    sz1      = floor(sz/2);
    [uu, vv] = meshgrid(-sz1(1):sz2(1), -sz1(2):sz2(2));
    if kernel_type == 6 % kernel for band-pass noise generation
      scale  = randi([1 4],1);
      dist   = (uu).^2 + (vv).^2;
      kernel = cos(sqrt(dist) / scale) .* fspecial('gaussian', [sz(1), sz(2)], 10);
    
    elseif kernel_type == 7 % kernel for line pattern noise generation
      scale  = randi([1 4],1);
      kernel = cos((uu + vv) / scale) .* fspecial('gaussian', [sz(1), sz(2)], 10);
    
    else % kernel for pink noise generation
      dist   = (uu).^2 + (vv).^2;
      spec3  = sqrt((sqrt(prod(sz))*1e-2)./(sqrt(dist) +  sqrt(prod(sz))*1e-2));
      kernel = fftshift(ifft2(ifftshift(spec3)));
    end
    kernel   = kernel / norm(kernel(:));
    
else
    
    kernel   = 1; % Dirac delta kernel for uncorrelated noise
    
end

kernel_power = min( max(rand/2,.05), .5); % to avoid insignificant or very strong noise
kernel = kernel_power*kernel; % convolutional kernel

%%% generating noise as circular convolution of standard white Gaussian noise with a kernel
if mod(size(kernel,1),2)==0
    noise = convn(padarray(normalizer(randn(SizeX)),(size(kernel)-1),'both','circular'),...
                  kernel(end:-1:1,end:-1:1),'valid');
    noise = noise(size(kernel,1)/2+1:end-size(kernel,1)/2+1,size(kernel,2)/2+1:end-size(kernel,2)/2+1);              
else
    noise = convn(padarray(normalizer(randn(SizeX)),(size(kernel)-1)/2,'both','circular'),...
                  kernel(end:-1:1,end:-1:1),'valid');
end
%%%%% equivalent to:
%%%% noise  = imfilter(normalizer(randn(SizeX)), kernel, 'circular');


%%% computing the FFT-PSD of the noise
noise_FFT_PSD = abs(fft2(kernel,SizeX(1),SizeX(2)).^2)*prod(SizeX);

end