%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point eigenvalues of pre-computed spiral wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%% Select system to solve
file_names.spiral_file = 'data/Bar_Eiswirth_spiral_r320_nx64_ny2400.mat';
file_names.out_file = 'Bar_Eiswirth_spectrum_eigs_r320_nx64_ny2400_wt';

%% Set up
load(file_names.spiral_file);
addpath ../utilities/

%% Standard linear operators
% file_names.jacobian_fcn ='jacobian_for_spectra_bar_eiswirth';
% [L1, L2] = ComputeLinearOperator_shortGrid(par,numPar);
% fh_jacobian = str2func(file_names.jacobian_fcn);
% [f,J] = fh_jacobian(U,L1,L2,par,numPar);

%% Weighted linear operators
par.w = 0.3; % weight used in eigs spectrum; set to zero for no weight
numVals = [50; 50; 2];
seeds = [-0.2-1i*0.45; -0.5-1i*0.25; -0.06-1i*0.54];
[L1, L2, L1r] = ComputeLinearOperator_shortGrid_with_Lr(par,numPar);
[~, J] = Bar_Eiswirth_weighted_operator(U, L1, L2, L1r, par, numPar);

%% Spectra
evals = [];
for k=1:length(numVals)
	[V, tmp, flag] = eigs(J,numVals(k),seeds(k));
	disp(flag)
	tmp = sort(diag(tmp),'ComparisonMethod','real');
	evals = [evals; tmp];
	plot(evals, '*');
end

save([file_names.out_file,'.mat'], 'evals')
evals = [real(evals), imag(evals)];
save([file_names.out_file,'.txt'], 'evals', '-ascii')
