%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve for a spiral wave in a given reaction-diffusion model
% Assumes you have an initial guess
% Assumes that r2 = 1 (radius scaled to be 1)
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify systems to solve and initial data
close all; clear all;

% Load files
file_names.spiral =   'data/Bar_Eiswirth_spiral_r40_nx64_ny300.mat'; % Initial data
file_names.out_name = 'Bar_Eiswirth_spiral_r40_nx64_ny300.mat';  % Out file 

file_names.solver_fcn = 'Bar_Eiswirth_2D_spiral_neumann';  % System to solve: F(U) = 0: can input Neumann or non-reflecting conditions here.
file_names.set_up_phase_condition = 'spiral_2D_phase_condition'; % Defines phase condition: modify for Neumann or non-reflecting conditions

%% Load data and set parameters
u0 = load(file_names.spiral);
U0 = u0.U;          % Initial guess
par = u0.par;       % Structure of model parameters
numPar = u0.numPar; % Structure of numerical parameters

% Set any parameter values (cannot be too different from initial condition - use continuation code for large changes)
par.Free = 'omega';             % Free variable (updated in solver step)
par.delta = 0.1;

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500);   % Options for fsolve 
addpath ../utilities/
fh_solver = str2func(file_names.solver_fcn);    % Set function handle for solver function

% Domain dependent setp: linear operators and phase condition
% Saves linear operators into the numPar struct
eval(file_names.set_up_phase_condition)

%% Solve for pattern
initial_sol = [U0; par.(par.Free)];
uout = fsolve(@(y) fh_solver(y,par,numPar,phase_cond),initial_sol,options); 
par.(par.Free) = uout(end);

% plot the pattern
plot_spiral(uout,par,numPar);

U = uout(1:par.numVars*numPar.nx*numPar.ny);
save(file_names.out_name ,'U','par','numPar');

