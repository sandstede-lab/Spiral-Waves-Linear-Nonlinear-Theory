%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically continue a spiral wave in a parameter value
% Assumes you have an initial guess
% Assumes that r2 = 1 (radius scaled to be 1)
% Works for model parameters and non-reflecting boundary conditions
% (depends on problem function)
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;

%% Select system to solve
% Modify for Rossler or Karma and Neumann or Non-reflectings
file_names.spiral =   'data/Barkley_spiral_r20_b0p05_a0p8.mat'; % Initial data
file_names.out_name = 'Barkley_spiral_r20_bm0p05_a0p8';  % Out file
file_names.problem = 'Barkley_2D_spiral_neumann';  % System to solve
file_names.set_up_phase_condition = 'spiral_2D_phase_condition'; % Defines phase condition and sets up linear operators

%% Set up
u0 = load(file_names.spiral);
U0 = u0.U;              % Initial solution
par = u0.par;           % Structure of model parameters
numPar = u0.numPar;     % Structure of numerical parameters
par.kappa = 0;

% Continuation parameter: picked from system parameter
contPar.final = -0.05;
contPar.ds = 2;				% continuation step size: this is for a secant continuation, so does not directly compare to step size of parameter. should always be positive
contPar.initial_ds = -0.01;		% initial step: sign determines direction of continuation
contPar.Name = 'b';				% Continuation parameter
contPar.Free = 'omega';			% Free variable
contPar.numContSteps = 200;		% Maximum number of continuation steps
contPar.plot_iter = 10;			% Plots continuation diagram and spiral every contPar.plot_iter steps

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...  % fsolve options
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500);

addpath ../utilities/
fh = str2func(file_names.problem);

%% Compute First two points
eval(file_names.set_up_phase_condition) % Set up the phase condition

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
free_vec = zeros(contPar.numContSteps,1);
cont_vec = zeros(contPar.numContSteps,1);
free_vec(1) = par.(contPar.Free);
cont_vec(1) = par.(contPar.Name);

%% Continuation
con_start = tic;

for k = 1:contPar.numContSteps
	disp(k)
	disp(['Continuation Parameter: ' num2str(par.(contPar.Name))])
	disp(['Free: ' num2str(par.(contPar.Free))])
	tmp_start = tic;
	
	% Predictor
	pred = sol1 + (sol1 - sol0)/norm(sol1 - sol0, 2) * contPar.ds;
	
	% Corrector - compute linear operator in the corrector
	uout = fsolve( @(sol) FixedPointSecantPredictorCorrector_general(sol,sol1,sol0,par,numPar,contPar,phase_cond,fh),pred,options);
	par.(contPar.Free) = uout(end-1);
	par.(contPar.Name) = uout(end);
	
	free_vec(k+1) = par.(contPar.Free);
	cont_vec(k+1) = par.(contPar.Name);
	
	% Prepare for next continuation step
	us = (phase_cond.pc).*uout(1:numPar.nx*numPar.ny);  % phase condition for the next round of iteration
	phase_cond.u_star = us(phase_cond.pc>0);
	
	us_th = (numPar.L1*uout(1:numPar.nx*numPar.ny));           % Will need to recompute linear operators if changing radius -- need to look into this (probably should be in the FixedPointSecantCorrector)
	phase_cond.u_star_th = us_th(phase_cond.pc > 0);
	
	sol0 = sol1;
	sol1 = uout;
	
	tmp_end = toc(tmp_start);
	disp(['fsolve Time: ' num2str(tmp_end)]);  % fsolve time
	
	if mod(k,contPar.plot_iter) == 0
		figure(4); plot(cont_vec(1:k+1),free_vec(1:k+1),'-o','LineWidth',2); % continuation diagram
		xlabel(contPar.Name,'FontSize',14); ylabel(contPar.Free,'FontSize',14); title('Continuation Parameter','FontSize',14); drawnow;
		% Spiral
		plot_spiral(uout,par,numPar);
		
	end
	
	% Stopping criteria: Not the greatest way to implement, but
	if (par.(contPar.Name)*sign(contPar.initial_ds)) >= (contPar.final*sign(contPar.initial_ds))
		disp(['Continuation parameter reached before specified number of iterations. Final parameter value: ' num2str(par.(contPar.Name))])
		break;
	end	
end

if k == contPar.numContSteps
	disp(['Iterations completed before continuation parameter reached. Current parameter value: ' num2str(par.(contPar.Name)) ' Desired value: ' num2str(contPar.final)])
else  % Do one final solve to set parameter value to desired result
	par.(contPar.Name) = contPar.final;
	uout = fsolve(@(y) fh(y,par,numPar,phase_cond),uout(1:end-1),options);
	par.(contPar.Free) = uout(end);
	cont_vec = cont_vec(1:k+1);
	free_vec = free_vec(1:k+1);
end

% Final continuation plot

figure(4); plot(cont_vec,free_vec,'-o','LineWidth',2);
xlabel(contPar.Name); ylabel(contPar.Free);
title('Continuation Parameter');
set(gca,'fontsize',16);
plot_spiral(uout,par,numPar);

% save final data
U = uout(1:par.numVars*numPar.nx*numPar.ny);
cont_data = [cont_vec, free_vec];

save([file_names.out_name,'.mat'] ,'U','par','numPar','cont_vec','free_vec','contPar');
save([file_names.out_name,'.txt'],'cont_data','-ascii');
save([file_names.out_name,'_start.txt'],'U0','-ascii');
save([file_names.out_name,'_end.txt'],'U','-ascii');

