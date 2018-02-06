%tvreg_demo1  Demo script for a TV tomography problem
%
% This script illustrates the use of the TV algorithm
% implemented in the functions tvreg_upn and tvreg_gpbb.
% The user can easily modify the script for
% other noise levels, regularization etc..
%
% The scripts loads a a clean threedimensional version of the 
% classical Shepp-Logan phantom image, and  
% obtains the observed data by multiplication with 
% a tomography matrix A and addition of noise.
% The algorithms tvreg_upn and tvreg_gpbb are used to 
% obtain the TV reconstructions of the phantom. An reference
% solution is calculated to obtain an estimate of the optimal
% objective.

clear
close all
clc

addpath lebdir tomobox

% Ensure we can reproduce results by fixing random numbers
randn('state',600)
rand('state',700)

%% Set parameters
rnl    = 0.01;       % Noise level
alpha  = 0.1;        % Regularization parameter
r1_max = 13;         % Halfwidth of object cube
N      = 2*r1_max+1; % Full width of object cube
N3     = N^3;        % Total number of variables
dims   = [N,N,N];    % Dimensions
u_max  = 15;         % Halfwidth of projection planes
U      = 2*u_max+1;  % Full width of projection planes
numProjections = 25; % Number of projections to reconstruct from

%% Construct tomography system matrix A and make spy plot
P              = getLebedevDirections(numProjections);
A              = buildSystemMatrix(r1_max,u_max,P);

figure
spy(A)
title('Sparsity pattern in A')

%% Construct true image and right hand side
X0     = phantom3d('Modified Shepp-Logan',N);       % Cube version
x0     = X0(:);                                     % Vectorized version
borig  = A*x0;                                      % Projections from true
Borig  = reshape(borig,U,U,numProjections);         % Projections in layers

%% Add noise
e = getNoise(rnl,borig);                % Gaussian white noise
b = borig+e;                            % Additive
B = reshape(b,U,U,numProjections);      % Noisy projections in layers

%% Display layers of true image along with noisefree and noisy projeections
figure
plotLayers(X0)
suptitle('True image')

figure
plotLayers(Borig)
suptitle('Noisefree projections')

figure
plotLayers(B)
suptitle('Noisy projections')

%% Parameters for the reconstruction algorithms
tau         = 1e-4*norm(x0,'inf');       % Huber smoothing parameter

% Specify nonnegativity constraints
constraint.type = 2;
constraint.c    = 0*ones(prod(dims),1);
constraint.d    = 1*ones(prod(dims),1);

% Options
opt.epsb_rel = 1e-6;
opt.k_max    = 10000;
opt.qs       = 1;
opt.K        = 2;
opt.verbose  = 1;
opt.beta     = 0.95;

% Options for reference solution
opt_ref.epsb_rel = 1e-8;
opt_ref.k_max    = 20000;
opt_ref.verbose  = 1;

%% Solve: Compute TV minimizer

% Reference solution
[x_ref fxk_ref hxk_ref gxk_ref fxkl_ref info_ref] = ...
    tvreg_upn(A,b,alpha,tau,dims,constraint,opt_ref);
fs = fxkl_ref(end);     % Final reference objective function value

% Solve using GPBB
tic
[xk_GPBB fxk_GPBB hxk_GPBB gx_kGPBB fxkl_GPBB info_GPBB] = ...
    tvreg_gpbb(A,b,alpha,tau,dims,constraint,opt);
tGPBB = toc

% Solve using UPN
tic
[xk_UPN fxk_UPN hxk_UPN gxk_UPN fxkl_UPN info_UPN] = ...
    tvreg_upn(A,b,alpha,tau,dims,constraint,opt);
tupn = toc

% Solve using UPN0
tic
opt.qs = 0;
[xk_UPNz fxk_UPNz hxk_UPNz gxk_UPNz fxkl_UPNz info_UPNz] = ...
    tvreg_upn(A,b,alpha,tau,dims,constraint,opt);
tupnz = toc

%% Plot convergence rates in terms of objective function values of the 
%% three methods, comparing to the final reference objective function value
figure
stairs(abs((fxkl_UPN-fs)/fs),'r')
hold on
stairs(abs((fxkl_UPNz-fs)/fs),'b')
stairs(abs((fxkl_GPBB-fs)/fs),'g')
set(gca,'yscale','log')
legend('UPN','UPN_0','GPBB')
xlabel('k')
ylabel('(f(x^k)-f^*)/f^*')
title('Convergence')

%% Display reconstructions
figure
plotLayers(reshape(xk_GPBB,dims))
suptitle('GPBB reconstruction')

figure
plotLayers(reshape(xk_UPN,dims))
suptitle('UPN reconstruction')

figure
plotLayers(reshape(xk_UPNz,dims))
suptitle('UPN_0 reconstruction')