%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point eigenvalues of pre-computed spiral wave
% Uses sparse methods (eigs) with seed points defined by 'seed_real' and
% 'seed_imag'
% Works for neumann or non-reflecting boundary conditions, depending on
% Jacobian function
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%% Select system to solve
% Modify for Karma or Rossler and non-reflecting or neumann boundary
% conditions
file_names.spiral_file = 'data/Rossler_spiral_R125_c2.mat';
%file_names.spiral_file = 'data/Rossler_spiral_R62_c2.mat';
%file_names.spiral_file = 'data/Rossler_spiral_R31_c2.mat';

file_names.out_file = 'Rossler_spectrum_R125_c2.txt';
%file_names.out_file = 'Rossler_spectrum_R62_c2.txt';
%file_names.out_file = 'Rossler_spectrum_R31_c2.txt';

file_names.jacobian_fcn = 'jacobian_for_spectra_rossler';  % Function to compute Jacobian: uses short grid. Function defines Neumann or Non-reflecting boundary conditions and Rossler or Karma

%% Set up
load(file_names.spiral_file);

% full point spectrum
numVals = [10;10;12;20;  30;  26;51;  30;  10;51];
seeds = [-0.35; -0.25; -0.15; -0.05
	-0.3+i*abs(par.omega)/2; ...
	-0.3+i*abs(par.omega); i*abs(par.omega); ...
	-0.3+i*abs(par.omega)*3/2; ...
	-0.3+i*abs(par.omega)*2; i*abs(par.omega)*2; ...
	];

% spectrum near period-doubling branch
%numVals = [30];
%seeds = [-0.3+i*abs(par.omega)/2];

addpath ../utilities/

[L1, L2] = ComputeLinearOperator_shortGrid(par,numPar); % This function does not use domain size of 1
fh_jacobian = str2func(file_names.jacobian_fcn);
[f,J] = fh_jacobian(U,L1,L2,par,numPar);

%% Spectra
eigenvals = [];
for k=1:length(numVals)
	tmp = eigs(J,numVals(k),seeds(k));
	tmp = sort(tmp,'ComparisonMethod','real');
	evals{k} = [real(tmp), imag(tmp)];
	eigenvals = [eigenvals; real(tmp), imag(tmp); NaN, NaN];
end

%% Plot
figure;
hold on;
for k=1:length(numVals)
	plot(evals{k}(:,1),evals{k}(:,2),'.','MarkerSize',18);
end
hold off;
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
box on;
set(gca,'fontsize',20,'linewidth',2);

save(file_names.out_file, 'eigenvals', '-ascii')
