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

u = load('rician_denoise');

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


figure(2);
imagesc(u);
axis image
axis off
colormap(gray(256));
title(sprintf('Denoised (PSNR %.2f dB)', ...
    -10*log10(mean((uexact(:) - u(:)).^2))));
