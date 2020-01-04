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
% Modify for Karma or Rossler and non-reflecting or neumann boundary conditions
file_names.spiral_file = 'data/Karma_spiral_R5_re1p2.mat';   % Precomputed spiral file
file_names.jacobian_fcn = 'jacobian_for_spectra_karma';  % Function to compute Jacobian: uses short grid. Function defines Neumann or Non-reflecting boundary conditions and Rossler or Karma
file_names.out_file = 'Karma_spectrum_point_R5_re1p2.txt';

%% Set up
load(file_names.spiral_file);

% Define the real and imaginary values for seed points
omega = abs(par.omega)
numVals = [5, 5, 10, 10, 10, 10];
seeds = [0.01, 0.01+i*omega, ...
	-10+i*omega/2, i*omega/2, ...
	-10+i*omega*3/2, i*omega*3/2];

addpath ../utilities/
[L1, L2] = ComputeLinearOperator_shortGrid(par,numPar); % This function does not use domain size of 1
fh_jacobian = str2func(file_names.jacobian_fcn);
[f,J] = fh_jacobian(U,L1,L2,par,numPar);

%% Spectra
eigenvals = [];
for k=1:length(numVals)
	[vecs{k}, tmp, flag] = eigs(J, numVals(k), seeds(k));
	disp(flag)
	tmp = sort(diag(tmp),'ComparisonMethod','real');
	evals{k} = [real(tmp), imag(tmp)];
	eigenvals = [eigenvals; real(tmp), imag(tmp); NaN, NaN];
end

%eigenvals = eig(full(J));

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
