%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve for a 2pi-periodic wave train in moving coordinate
% frame
% Assumes you have an initial guess
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; 

%% System to solve
% Modify for Karma or Rossler
file_names.wave_train = 'data/Barkley_wave_train_b0p05_a0p8.mat';				% Initial guess
par_load = load('data/Barkley_spiral_r20_bcont_a0p8.txt');
file_names.out_file   = 'Barkley_spiral_dispersion_r20_bcont_a0p8.txt';		% Name of output saved data
file_names.wt_solver = 'Barkley_1D';      % Function to solve the 1D wave train

%% Set up
b = par_load(:,1);
omega = par_load(:,2);

u0 = load(file_names.wave_train);
U0 = u0.U;
numPar = u0.numPar;
par = u0.par;

% Set the free parameter
par.Free = 'kappa';
V0 = [U0; par.(par.Free)];

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',... % fsolve options
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',1500);      
addpath ../utilities/
fh = str2func(file_names.wt_solver);

%% Do an fsolve to find solution

[L1, L2] = ComputeLinearOperator_1D(par,numPar);
for k=1:length(b)
	par.b = b(k);
	par.omega = omega(k);
	phase_cond.u_old = V0(1:numPar.nx);

	V0 = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),V0,options);  
	par.(par.Free) = V0(end);
	par_load(k,3) = V0(end);
end

save(file_names.out_file,'par_load','-ascii');
