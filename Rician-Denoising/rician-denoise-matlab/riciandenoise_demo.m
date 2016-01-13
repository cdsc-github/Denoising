%%% Demo script for riciandenoise.m %%%


% Get a clean input image
uexact = imread('knee-mri.png');
uexact = mean(double(uexact)/255,3);
% Simulate Rician noise
sigma = 0.05;
f = ricianrnd(uexact,sigma);

% Parameters for riciandenoise
% Smaller lambda implies stronger denoising
lambda = 0.065;
% Parameter Tol is the stopping tolerance
Tol = 2e-3;
% Parameter FastApprox decides how to approximate I1(x)/I0(x): 
FastApprox = true;   % Fast and moderately accurate
% FastApprox = false;  % Slow but very accurate


% Plot input
figure(1);
imagesc(f);
axis image
axis off
colormap(gray(256));
title(sprintf('Noisy input image (PSNR %.2f dB)',...
    -10*log10(mean((uexact(:) - f(:)).^2))));
drawnow;
shg;

% Perform the denoising
StartTime = clock;
u = riciandenoise(f,sigma,lambda,Tol,FastApprox);
%u = riciandenoisemx(f,sigma,lambda,Tol);
StopTime = clock;

% Plot output
figure(2);
imagesc(u);
axis image
axis off
colormap(gray(256));
title(sprintf('Denoised (PSNR %.2f dB, CPU time %.2f s)', ...
    -10*log10(mean((uexact(:) - u(:)).^2)), ...
    etime(StopTime,StartTime)));
