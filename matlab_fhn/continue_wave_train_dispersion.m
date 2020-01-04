%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve for a 2pi-periodic wave train in moving coordinate frame
% Compute kappa and d_perp given deltra and omega from spiral dispersion relation
% Stephanie Dodson & Bjorn Sandstede
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; 

%% System to solve
% Modify for Karma or Rossler
file_names.wave_train = 'data/FHN_wave_train_1024_delta1p1.mat';  % Initial guess
par_load = load('data/FHN_spiral_continuation.txt');  % spiral dispersion relation
file_names.out_file   = 'FHN_spiral_continuation.txt';	% Name of output saved data
file_names.wt_solver = 'FHN_1D';   % Function to solve the 1D wave train

%% Set up
delta = par_load(:,1);
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
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',50);      
addpath ../utilities/
fh = str2func(file_names.wt_solver);

%% Do an fsolve to find solution

[L1, L2] = ComputeLinearOperator_1D(par,numPar);
for k=1:length(delta)
	par.delta = delta(k);
	par.omega = omega(k);
	phase_cond.u_old = V0(1:numPar.nx);

	V0 = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),V0,options);  
	par.(par.Free) = V0(end);
	par_load(k,3) = V0(end);

	U0 = V0(1:end-1);
 	J = FHN_1D_adjoint(U0,L1,L2,par,numPar);
 	[uadj, tmp] = eigs(J', 1, 0);
 	u_prime = [L1*U0(1:end/2); L1*U0(end/2+1:end)];
 	Du_prime = [u_prime(1:end/2); par.delta*u_prime(end/2+1:end)];
 	par_load(k,4) = (uadj'*Du_prime)/(uadj'*u_prime);
end

save(file_names.out_file,'par_load','-ascii');
