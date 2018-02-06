%tvreg_demo2  Demo script for a TV deblurring problem
%
% This script illustrates the use of the TV algorithm
% implemented in the functions tvreg_upn and tvreg_gpbb.
% The user can easily modify the script for
% other noise levels, regularization etc..
%
% The scripts loads a a clean image and  
% obtains the observed blurred image by convolution with 
% a point spread function and adding noise.
% The algorithms tvreg_upn and tvreg_gpbb are used to 
% obtain the TV reconstructions of the original image. A reference
% solution is calculated to obtain an estimate of the optimal
% objective.

clear
close all
clc

addpath tomobox

% Ensure we can reproduce results by fixing random numbers
randn('state',600)
rand('state',700)

%% Set parameters
rnl    = 0.01;      % Noise level
alpha  = 1.0;       % Regularization parameter

%% Set up point spread function (PSF) which describes the blurring
ks          = 7; 
ksigma      = 1.0;
c           = 4;
w           = exp(-((1:ks)-c).^2/(2*ksigma))';
PSF         = w*w';
PSF         = PSF/norm(PSF,'fro');
PSFs.PSF    = PSF;
PSFs.center = [c,c];

%% Load image and extract region of interest (ROI)
X0full = double(imread('Pirate.tif'));
m      = 128;
n      = 128;
offset = 10;
X0     = X0full(offset:m+offset-1,offset:n+offset-1);

%% Obtain blurred image b by convolving X0 with PSF. Extract ROI. Add noise
Bfull = conv2(X0full,PSFs.PSF,'same');
Borig = Bfull(offset:m+offset-1,offset:n+offset-1);
borig = Borig(:);
e     = getNoise(rnl,borig);
b     = borig+e;
B     = reshape(b,[m,n]);

%% Parameters for the reconstruction algorithms
tau  = 1e-4*norm(X0(:),'inf');
dims = [m,n];

% Specify box constraints [0,255]
constraint.type = 2;
constraint.c    = 0*ones(prod(dims),1);
constraint.d    = 255*ones(prod(dims),1);

% Options
opt.epsb_rel = 1e-4;
opt.k_max    = 10000;
opt.x0       = b;
opt.verbose  = 1;

% Options for reference solution
opt_ref.epsb_rel = 1e-6;
opt_ref.k_max    = 10000;
opt_ref.x0       = b;
opt_ref.verbose  = 1;

%% Solve: Compute TV minimizer

% Reference solution
[x_ref fxk_ref hxk_ref gxk_ref fxkl_ref] = ...
    tvreg_upn(PSFs,b,alpha,tau,dims,constraint,opt_ref);
fs = fxkl_ref(end);     % Final reference objective function value

% Solve using GPBB
tic
[xk_GPBB fxk_GPBB hxk_GPBB gx_kGPBB fxkl_GPBB] = ...
    tvreg_gpbb(PSFs,b,alpha,tau,dims,constraint,opt);
GPBB = toc;

% Solve using UPN
tic
opt.qs=1;
[xk_UPN fxk_UPN hxk_UPN gxk_UPN fxkl_UPN] = ...
    tvreg_upn(PSFs,b,alpha,tau,dims,constraint,opt);
tUPN = toc

% Solve using UPN0
tic
opt.qs=0;
[xk_UPNz fxk_UPNz hxk_UPNz gxk_UPNz fxkl_UPNz] = ...
    tvreg_upn(PSFs,b,alpha,tau,dims,constraint,opt);
tUPN0 = toc

%% Plot convergence rates in terms of objective function values of the 
%% three methods, comparing to the final reference objective function value
figure
clf
stairs(abs((fxkl_UPN-fs)/fs),'r')
hold on
stairs(abs((fxkl_UPNz-fs)/fs),'b')
stairs(abs((fxkl_GPBB-fs)/fs),'g')
set(gca,'yscale','log')
legend('UPN','UPN_0','GPBB')
xlabel('k')
ylabel('(f(x^k)-f^*)/f^*')

%% Display reconstructions
figure
colormap gray

subplot(2,3,1)
imagesc(X0);
axis image
title('Original')

subplot(2,3,2)
imagesc(Borig);
axis image
title('Blurred')

subplot(2,3,3)
imagesc(B);
axis image
title('Blurred and noisy')

subplot(2,3,4)
imagesc(reshape(xk_GPBB,m,n));
axis image
title('GPBB reconstruction')

subplot(2,3,5)
imagesc(reshape(xk_UPN,m,n));
axis image
title('UPN reconstruction')

subplot(2,3,6)
imagesc(reshape(xk_UPN,m,n));
axis image
title('UPN_0 reconstruction')