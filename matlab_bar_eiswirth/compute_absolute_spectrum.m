%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes absolute spectrum of spiral wave
% Assumes an initial starting point for the absolute spectrum
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
file_names.starting_point = 'Bar_Eiswirth_spatial_spectrum_nx256.mat';
file_names.out_name = 'Bar_Eiswirth_absolute_spectrum_nx256.txt';
file_names.cont_func = 'absolute_spectra_bar_eiswirth';
file_names.spatial_evals_fcn = 'spatial_evals_fcn_bar_eiswirth';

%%
soln = load(file_names.starting_point);
u_infty = soln.u_infty;                  % Asymptotic wave trains
numPar = soln.numPar;
par = soln.par;

U0 = [soln.uout(1:numPar.nx); soln.uout(2*numPar.nx+1:3*numPar.nx)]; % Eigenvector corresponding to par.nu
W0 = zeros(size(U0));              % Second eigenvectors (of second spatial eigenvalue par.nu + 1i*par.beta) - random guesses, which will be solved for

contPar.ds = c_direction*0.001;           % Continuation step size: sign determines direction of continuation
contPar.Name = 'beta';       % Continuation parameter: vertical separation between spatial eigenvalues: Im(nu_1) - Im(nu_1)
contPar.free1 = 'nu';        % Free parameter 1
contPar.free2 = 'lambda';    % Free parameter 2
contPar.numContSteps = 5000;
contPar.sp_evals_iter = 5000;   % Computes spatial eigenvalues and updates plots ever contPar.sp_evals_iter iterations

lambda_start = par.lambda;
beta_start = par.beta;

% Additional variables
options = optimset('Display','none','Jacobian','on', 'DerivativeCheck','off',...
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
cont_vec = zeros(contPar.numContSteps+1,1);		% beta
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

	if mod(k,contPar.sp_evals_iter) == 0
		[A,V,vals] = fh_spatial(L1,L2,u_infty,par,numPar);
		[vals, idx] = sort(vals,'ComparisonMethod','real');
		if (abs(vals(end/2)-par.(contPar.free1))>0.01) & (abs(vals(end/2+1)-par.(contPar.free1))>0.01)
			disp('Morse index changed');
			vals(end/2)
			vals(end/2+1)
			par.(contPar.free1)
			free1_vec = free1_vec(1:k);
			free2_vec = free2_vec(1:k);
			break;
		end
	end

	if real(par.beta)*real(beta_start)<0
		disp('Reached branch point');
		free1_vec = free1_vec(1:k);
		free2_vec = free2_vec(1:k);
		break
	end
end

if c_flip==1
	evals = [real(free2_vec(end:-1:1)), imag(free2_vec(end:-1:1))];
else
	evals = [real(free2_vec), imag(free2_vec)];
end

if c_append==1
	save(file_names.out_name, 'evals', '-ascii', '-append');
else
	save(file_names.out_name, 'evals', '-ascii');
end

