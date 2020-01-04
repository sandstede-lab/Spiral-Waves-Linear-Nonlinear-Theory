%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically continue a spiral wave in a parameter value
%     and computes leading spectrum alongside
% Assumes you have an initial guess
% Assumes that r2 = 1 (radius scaled to be 1)
% Works for model parameters and non-reflecting boundary conditions
%     (depends on problem function)
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

%% select system to solve
file_names.spiral =   'data/Barkley_spiral_outward.mat'; % Initial data
file_names.spectra = 'Barkley_spectrum_outward.txt';

file_names.problem = 'Barkley_2D_spiral_neumann';  % System to solve
file_names.set_up_phase_condition = 'spiral_2D_phase_condition'; % Defines phase condition and sets up linear operators
file_names.jacobian_fcn = 'jacobian_for_spectra_barkley';

%% set up continuation
u0 = load(file_names.spiral);
U0 = u0.U;              % Initial solution
par = u0.par;           % Structure of model parameters
numPar = u0.numPar;     % Structure of numerical parameters
par.kappa = 0;

% continuation parameter: picked from system parameter

%% inward
%contPar.final = 0.58;
%par.a = 0.535;
%seed = [0, 1i*abs(par.omega), -0.1716+1i*0.9020];

%% outward
contPar.final = 0.76;
par.a = 0.72;
seed = [0, 1i*abs(par.omega), 0.37+1i*1.94];

contPar.ds = 0.5;				% continuation step size: this is for a secant continuation, so does not directly compare to step size of parameter. should always be positive
contPar.initial_ds = 0.005;		% initial step: sign determines direction of continuation
contPar.Name = 'a';				% Continuation parameter
contPar.Free = 'omega';			% Free variable
contPar.numContSteps = 200;		% Maximum number of continuation steps
contPar.plot_iter = 200;			% Plots continuation diagram and spiral every contPar.plot_iter steps

% additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...  % fsolve options
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500);

addpath ../utilities/

% setup continuation
fh = str2func(file_names.problem);
eval(file_names.set_up_phase_condition) % Set up the phase condition

% set up spectra
[L1, L2] = ComputeLinearOperator_shortGrid(par,numPar); % This function does not use domain size of 1
fh_jacobian = str2func(file_names.jacobian_fcn);

%% Compute First two points
initial_sol = [U0; par.(contPar.Free)];
uout = fsolve(@(y) fh(y,par,numPar,phase_cond),initial_sol,options);

par.(contPar.Free) = uout(end);
sol0 = [uout; par.(contPar.Name)];   % First point

% Phase condition for next step
us = (phase_cond.pc).*uout(1:numPar.nx* numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);
us_th = (numPar.L1*uout(1:numPar.nx* numPar.ny)); % d/d(theta) of phase condition
phase_cond.u_star_th = us_th(phase_cond.pc > 0);

% Compute the second point
par.(contPar.Name) = par.(contPar.Name) + contPar.initial_ds;   % Add initial step to continuation parameter

% Solve system to compute second point: can be slow - depends on contPar.ds
uout = fsolve(@(y) fh(y,par,numPar,phase_cond),uout,options);
par.(contPar.Free) = uout(end);
sol1 = [uout; par.(contPar.Name)]; % Second point

us = (phase_cond.pc).*uout(1:numPar.nx*numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);
us_th = (numPar.L1*uout(1:numPar.nx*numPar.ny));
phase_cond.u_star_th = us_th(phase_cond.pc > 0);

%% Continuation diagram
results = [];

%% Continuation
con_start = tic;

for k = 1:contPar.numContSteps
	disp(k)
	disp(['Continuation Parameter: ' num2str(par.(contPar.Name))])
	disp(['Free: ' num2str(par.(contPar.Free))])
	tmp_start = tic;
	
	% Predictor
	pred = sol1 + (sol1 - sol0)/norm(sol1 - sol0, 2) * contPar.ds;
	
	% Corrector
	uout = fsolve( @(sol) FixedPointSecantPredictorCorrector_general(sol,sol1,sol0,par,numPar,contPar,phase_cond,fh),pred,options);
	par.(contPar.Free) = uout(end-1);
	par.(contPar.Name) = uout(end);

	% Spectra
	[f,J] = fh_jacobian(uout(1:end-1),L1,L2,par,numPar);

	lambda = zeros(1,3);
	for k=1:length(seed)
		lambda(k) = eigs(J,1,seed(k));
	end
	results = [results; par.a, real(lambda), imag(lambda)];
	seed = lambda;

	% Prepare for next continuation step
	us = (phase_cond.pc).*uout(1:numPar.nx*numPar.ny);  % phase condition for the next round of iteration
	phase_cond.u_star = us(phase_cond.pc>0);
	
	us_th = (numPar.L1*uout(1:numPar.nx*numPar.ny));
	phase_cond.u_star_th = us_th(phase_cond.pc > 0);
	
	sol0 = sol1;
	sol1 = uout;
	
	tmp_end = toc(tmp_start);
	disp(['fsolve Time: ' num2str(tmp_end)]);  % fsolve time
	
	if (par.(contPar.Name)*sign(contPar.initial_ds)) >= (contPar.final*sign(contPar.initial_ds))
		disp(['Continuation parameter reached before specified number of iterations. Final parameter value: ' num2str(par.(contPar.Name))])
		break;
	end	
end

U = uout(1:par.numVars*numPar.nx*numPar.ny);
save(file_names.spectra,'results','-ascii');
	