%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point eigenvalues of adjoint of pre-computed spiral wave
% Uses sparse methods (eigs)
% Works for neumann or non-reflecting boundary conditions, depending on
% Jacobian function
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%% Select system to solve
file_names.spiral_file = 'data/Barkley_spiral_r20_bm0p05_a0p8.mat';
file_names.out_file = 'Barkley_adjoint_r20_bm0p05_a0p8.txt';

%% Set up system
file_names.jacobian_fcn = 'jacobian_for_spectra_barkley';
load(file_names.spiral_file);
addpath ../utilities/

[L1, L2] = ComputeLinearOperator_shortGrid(par,numPar);
fh_jacobian = str2func(file_names.jacobian_fcn);
[f,J] = fh_jacobian(U,L1,L2,par,numPar);
J = transpose(J);

%% Compute eigenvalues
seed = 0;

[V,D] = eigs(J,1,seed);
D

nx = numPar.nx;
ny = numPar.ny;
m = 1+nx*(ny-1);

v = abs(V);
V1 = [ones(nx,1)*v(1);     v(2:m)];
V2 = [ones(nx,1)*v(1+m);   v(2+m:2*m)];

S = [U, [V1; V2]];
save(file_names.out_file, 'S', '-ascii')

