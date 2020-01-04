%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point eigenvalues of adjoint of pre-computed spiral wave
% Uses sparse methods (eigs)
% Works for neumann or non-reflecting boundary conditions, depending on
% Jacobian function
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%% Select system to solve
file_names.spiral_file = 'data/Karma_spiral_R5_re1p2.mat';
file_names.out_file = 'Karma_adjoint_R5_re1p2.txt';

%% Set up system
file_names.jacobian_fcn = 'jacobian_for_spectra_karma';
load(file_names.spiral_file);
addpath ../utilities/

[L1, L2] = ComputeLinearOperator_shortGrid(par,numPar);
fh_jacobian = str2func(file_names.jacobian_fcn);
[f,J] = fh_jacobian(U,L1,L2,par,numPar);

nx = numPar.nx;
ny = numPar.ny;
m = 1+nx*(ny-1);

%% Compute adjoint eigenfunction at zero
Jadj = transpose(J);
seed = 0;
[V,D] = eigs(Jadj,1,seed);
D

v = abs(V);
V1 = [ones(nx,1)*v(1);     v(2:m)];
V2 = [ones(nx,1)*v(1+m);   v(2+m:2*m)];
S = [U, [V1; V2]];

%% Compute unstable eigenfunction
seed = 6.3195743e-01 + 1i*8.0722642e+01;
[V,D] = eigs(J,1,seed);
D

v = abs(V);
V1 = [ones(nx,1)*v(1);     v(2:m)];
V2 = [ones(nx,1)*v(1+m);   v(2+m:2*m)];
S = [S, [V1; V2]];

%% save solution
save(file_names.out_file, 'S', '-ascii')
