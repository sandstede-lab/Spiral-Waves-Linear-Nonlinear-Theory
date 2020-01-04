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
file_names.spiral =   'data/FHN_spiral_solved_nx128_ny400_r2_20.mat';
file_names.out_name = 'FHN_spiral_continuation.txt';

file_names.problem = 'FHN_2D_spiral_neumann';
file_names.set_up_phase_condition = 'spiral_2D_phase_condition';

%% Set up
u0 = load(file_names.spiral);
U0 = u0.U;              % Initial solution
par = u0.par;           % Structure of model parameters
numPar = u0.numPar;     % Structure of numerical parameters
%disp(par);
%disp(numPar);

% Continuation parameter: picked from system parameter
contPar.final = 1.16;
contPar.ds = 1;
contPar.initial_ds = 0.001;		% initial step: sign determines direction of continuation
contPar.Name = 'delta';			% Continuation parameter
contPar.Free = 'omega';			% Free variable
contPar.numContSteps = 1000;	% Maximum number of continuation steps
contPar.plot_iter = 10;			% Plots continuation diagram and spiral every contPar.plot_iter steps

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...  % fsolve options
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',50);

addpath ../utilities/
fh = str2func(file_names.problem);

%% Compute First two points
eval(file_names.set_up_phase_condition) % Set up the phase condition

initial_sol = [U0; par.(contPar.Free)];
uout = fsolve(@(y) fh(y,par,numPar,phase_cond),initial_sol,options);

par.(contPar.Free) = uout(end);
sol0 = [uout; par.(contPar.Name)];   % First point

cont = [par.(contPar.Name), par.(contPar.Free)];

% Phase condition for next step
us = (phase_cond.pc).*uout(1:numPar.nx* numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);
us_th = (numPar.L1*uout(1:numPar.nx* numPar.ny)); % d/d(theta) of phase condition
phase_cond.u_star_th = us_th(phase_cond.pc > 0);

% Compute the second point
par.(contPar.Name) = par.(contPar.Name)+contPar.initial_ds;   % Add initial step to continuation parameter

% Solve system to compute second point: can be slow - depends on contPar.ds
uout = fsolve(@(y) fh(y,par,numPar,phase_cond),uout,options);
par.(contPar.Free) = uout(end);
sol1 = [uout; par.(contPar.Name)]; % Second point
cont = [cont; par.(contPar.Name), par.(contPar.Free)];

us = (phase_cond.pc).*uout(1:numPar.nx*numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);
us_th = (numPar.L1*uout(1:numPar.nx*numPar.ny));
phase_cond.u_star_th = us_th(phase_cond.pc > 0);

%% Continuation diagram


%% Continuation
con_start = tic;

for k = 1:contPar.numContSteps
	disp(k)
	disp(['Continuation Parameter: ' num2str(par.(contPar.Name))])
	disp(['Free: ' num2str(par.(contPar.Free))])
	tmp_start = tic;
	
	% Predictor
	pred = sol1 + (sol1-sol0)/norm(sol1-sol0,2)*contPar.ds;
	
	% Corrector - compute linear operator in the corrector
	[uout,~,exitflag] = fsolve( @(sol) FixedPointSecantPredictorCorrector_general(sol,sol1,sol0,par,numPar,contPar,phase_cond,fh),pred,options);
	if exitflag<=0
		break;
	end
	
	par.(contPar.Free) = uout(end-1);
	par.(contPar.Name) = uout(end);
	cont = [cont; par.(contPar.Name), par.(contPar.Free)];
	
	% Prepare for next continuation step
	us = (phase_cond.pc).*uout(1:numPar.nx*numPar.ny);  % phase condition for the next round of iteration
	phase_cond.u_star = us(phase_cond.pc>0);
	
	us_th = (numPar.L1*uout(1:numPar.nx*numPar.ny));
	phase_cond.u_star_th = us_th(phase_cond.pc > 0);
	
	sol0 = sol1;
	sol1 = uout;
	
	tmp_end = toc(tmp_start);
	disp(['fsolve Time: ' num2str(tmp_end)]);  % fsolve time
	
	if mod(k,contPar.plot_iter) == 0
		% Continuation diagram
		figure(2);
		plot(cont(1:k+1,1),cont(1:k+1,2),'-o','LineWidth',2); % continuation diagram
		xlabel(contPar.Name,'FontSize',14); ylabel(contPar.Free,'FontSize',14); title('Continuation Parameter','FontSize',14);
		drawnow;
		% Spiral
		plot_spiral(uout,par,numPar);
		% Save
		U = uout(1:par.numVars*numPar.nx*numPar.ny);
		save(file_names.out_name,'U','par','numPar','cont');
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
	cont = [cont; par.(contPar.Name), par.(contPar.Free)];
end

% Final continuation plot
figure(2); plot(cont(:,1),cont(:,2),'-o','LineWidth',2);
xlabel(contPar.Name); ylabel(contPar.Free);
title('Continuation Parameter');
set(gca,'fontsize',16);
plot_spiral(uout,par,numPar);

% save final data
save(file_names.out_name, 'cont', '-ascii');
