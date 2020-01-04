%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically continue absolute spectrum of spiral wave
% Assumes you have an initial starting point for the absolute spectrum
% Uses the 2pi-periodic wave train of the spiral wave
% Continues along spectral curve and computes all spatial eigenvalues to
% ensure overlapping spatial eigenvalue criteria for absolute spectrum is met
% Simple continuation
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

c_direction = -1;
c_flip = 0;
c_append = 0;

%% Select system to solve
% Modify for Karma or Rossler
file_names.starting_point = 'data/Karma_start_abs_up_R5_re1p2.mat';
file_names.out_name = 'Karma_spectrum_absolute_R5_re1p2.txt';

file_names.cont_func = 'absolute_spectra_karma';
file_names.spatial_evals_fcn = 'spatial_evals_fcn_karma';

%%
soln = load(file_names.starting_point);
u_infty = soln.u_infty;                  % Asymptotic wave trains
numPar = soln.numPar;
par = soln.par;

U0 = soln.uout(1:par.numVars*numPar.nx); % Eigenvector corresponding to par.nu
W0 = zeros(size(U0));              % Second eigenvectors (of second spatial eigenvalue par.nu + 1i*par.beta) - random guesses, which will be solved for

contPar.ds = c_direction()*0.1;           % Continuation step size: sign determines direction of continuation
contPar.Name = 'beta';       % Continuation parameter: vertical separation between spatial eigenvalues: Im(nu_1) - Im(nu_1)
contPar.free1 = 'nu';        % Free parameter 1
contPar.free2 = 'lambda';    % Free parameter 2
contPar.numContSteps = 200;
contPar.sp_evals_iter = 5;   % Computes spatial eigenvalues and updates plots ever contPar.sp_evals_iter iterations

lambda_start = par.lambda;
beta_start = par.beta;

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',5000);

addpath ../utilities/
fh = str2func(file_names.cont_func);
fh_spatial = str2func(file_names.spatial_evals_fcn);

%% Initial data
% Compute linear operator
[L1, L2] = ComputeLinearOperator_1D(par,numPar); % 1st and 2nd derivative matrices: assumes 2pi-periodic domain

% Initial condition is the loaded data
initial_sol = [U0; W0; par.(contPar.free1); par.(contPar.free2)]; % Eigenvectors
phase_cond.u_old = U0;
phase_cond.w_old = W0;

% Do an fsolve to make sure eigenvector is right
uout = fsolve(@(y) fh(y,u_infty,L1,L2,par,numPar,phase_cond),initial_sol,options);
par.(contPar.free1) = uout(end-1);
par.(contPar.free2) = uout(end);

phase_cond.u_old = uout(1:par.numVars*numPar.nx);
phase_cond.w_old = uout(par.numVars*numPar.nx+1:2*par.numVars*numPar.nx);

%% Continuation
free1_vec = zeros(contPar.numContSteps+1,1);	% spatial eigenvalues
cont_vec = zeros(contPar.numContSteps+1,1);		% gamma
free2_vec = zeros(contPar.numContSteps+1,1);	% temporal eigenvalues
diff_vals = zeros(contPar.numContSteps,1);

free1_vec(1) = par.(contPar.free1);
free2_vec(1) = par.(contPar.free2);
cont_vec(1) = par.(contPar.Name);

for k = 1:contPar.numContSteps
	par.(contPar.Name) = par.(contPar.Name) + contPar.ds; % Update continuation parameter: simple continuation
	disp([contPar.Name ': ' num2str(par.(contPar.Name))])
	
	uout = fsolve(@(y) fh(y,u_infty,L1,L2,par,numPar,phase_cond),uout,options);
	par.(contPar.free1) = uout(end-1);
	par.(contPar.free2) = uout(end);
	phase_cond.u_old = uout(1:par.numVars*numPar.nx);
	phase_cond.w_old = uout(par.numVars*numPar.nx+1:2*par.numVars*numPar.nx);
	
	free1_vec(k+1) = par.(contPar.free1);
	free2_vec(k+1) = par.(contPar.free2);
	cont_vec(k+1) = par.(contPar.Name);
	
	% Additional stopping criteria
%	if real(par.lambda)<lambda_min
%		disp('Reached minimal lambda value');
%		free1_vec = free1_vec(1:k);
%		free2_vec = free2_vec(1:k);
%		break
%	end
%	if real(par.lambda-lambda_triple)*real(lambda_start-lambda_triple)<0
%		disp('Reached triple point');
%		free1_vec = free1_vec(1:k);
%		free2_vec = free2_vec(1:k);
%		break
%	end
	if real(par.beta)*real(beta_start)<0
		disp('Reached branch point');
		free1_vec = free1_vec(1:k);
		free2_vec = free2_vec(1:k);
		break
	end
end

if c_flip()==1
	evals = [real(free2_vec(end:-1:1)), imag(free2_vec(end:-1:1))];
else
	evals = [real(free2_vec), imag(free2_vec)];
end

evals = [evals; NaN, NaN; evals(:,1), par.omega-evals(:,2); NaN, NaN];
evals = [evals; evals + [0, par.omega]];

if c_append()==1
	save(file_names.out_name, 'evals', '-ascii', '-append');
else
	save(file_names.out_name, 'evals', '-ascii');
end

