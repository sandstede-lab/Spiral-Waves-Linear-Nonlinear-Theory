% Compute floquet and transverse spectrum using continuation methods outlined in
% Rademacher, Scheel, and Sandstede 2007
% Modify continuation (number of steps, branch, etc) through contPar struc

close all; clear all;

%% Load files
file_names.wt = 'data/FHN_spectrum_essential_start_delta1p16.mat';

file_names.out_name = 'FHN_spectrum_essential_delta1p16.txt';
%file_names.out_name = 'FHN_spectrum_transverse_delta1p16.txt';

%% Set up
file_names.ess_spec_solver = 'FHN_essSpec_fcn';    % Solve essential spectrum function: spiral dispersion relation

wt = load(file_names.wt);
par =    wt.par;
numPar = wt.numPar;
U0 = wt.ueig;
u_infty = wt.u_infty;         % Asymptotic wave train
par.lambda = 0;
par.nu = 0;

%% switch
par.flag = 1;	% =0 for transverse; =1 for floquet spectrum

% Continuation parameters for essential spectrum continuation -- This is a SIMPLE
% continuation!
contPar.ds   = 0.01*1i;
contPar.Name = 'nu';
contPar.free = 'lambda';
contPar.numContSteps = 100;
contPar.full_branch = 1;

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...  % Options for fsolve
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500);
addpath ../utilities/
ess_spec_solver_fh = str2func(file_names.ess_spec_solver);

%% Find initial point for continuation

% Compute linear operator
[L1, L2] = ComputeLinearOperator_1D(par,numPar); % 1st and 2nd derivative matrices in 1D

init_lambda = par.lambda;
init_nu = par.nu;
initial_sol = [U0; par.(contPar.free)]; % Initial eigenvectors
phase_cond.u_old = U0;

% Make sure eigenvector is correct
uout = fsolve(@(y) ess_spec_solver_fh(y,u_infty,L1,L2,par,numPar,phase_cond),initial_sol,options);
par.(contPar.free) = uout(end);
phase_cond.u_old = uout(1:par.numVars*numPar.nx);
initial_sol = uout;

%% Continuation
free_vec = zeros(contPar.numContSteps+1,1);
cont_vec = zeros(contPar.numContSteps+1,1);

free_vec(1) = par.lambda;
cont_vec(1) = par.nu;

for k = 1:contPar.numContSteps
	par.(contPar.Name) = par.(contPar.Name) + contPar.ds;   % Add step to continuation parameter
	
	uout = fsolve(@(y) ess_spec_solver_fh(y,u_infty,L1,L2,par,numPar,phase_cond),uout,options);
	par.(contPar.free) = uout(end);
	phase_cond.u_old = uout(1:par.numVars*numPar.nx);
	
	free_vec(k+1)   = par.(contPar.free);
	cont_vec(k+1)   = par.(contPar.Name);
end

%% Continuation in other direction to get full branch
if contPar.full_branch == 1
	free_vec2 = zeros(contPar.numContSteps,1);
	cont_vec2 = zeros(contPar.numContSteps,1);
	
	par.lambda = init_lambda;  % Reset to start at origin again
	par.nu = init_nu;
	uout = initial_sol;
	phase_cond.u_old = U0;
	
	contPar.ds = -1*(contPar.ds);
	
	for k = 1:contPar.numContSteps
		par.(contPar.Name) = par.(contPar.Name) + contPar.ds;   % Add step to continuation parameter
		
		uout = fsolve(@(y) ess_spec_solver_fh(y,u_infty,L1,L2,par,numPar,phase_cond),uout,options);
		par.(contPar.free) = uout(end);
		phase_cond.u_old = uout(1:par.numVars*numPar.nx);
		
		free_vec2(k)   = par.(contPar.free);
		cont_vec2(k)   = par.(contPar.Name);
	end
	
	ess_spec = [flip(free_vec); free_vec2];  % form full essential spectrum branch
	nu_vec   = [flip(cont_vec); cont_vec2];
else
	ess_spec = free_vec;
	nu_vec   = cont_vec;	
end

%% Plot and save data
figure(1);   % Plot the spectrum branch
subplot(1,2,1); hold on;
plot(real(ess_spec), imag(ess_spec), 'b*-', 'linewidth', 3);
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
title(['Essential Spectrum']);
plot([0,0], ylim,'color','k','linewidth',2);
plot(xlim, [0,0],'color','k','linewidth',2);
subplot(1,2,2)
plot(imag(nu_vec), real(ess_spec), 'b*-', 'linewidth', 3); hold on;
xlabel('Im(\nu)');  ylabel('Re(\lambda)');
title(['Essential Spectrum']);
plot([0,0], ylim,'color','k','linewidth',2);
plot(xlim, [0,0],'color','k','linewidth',2);

results = [real(ess_spec), imag(ess_spec), real(nu_vec), imag(nu_vec)];
save(file_names.out_name,'results','-ascii');
