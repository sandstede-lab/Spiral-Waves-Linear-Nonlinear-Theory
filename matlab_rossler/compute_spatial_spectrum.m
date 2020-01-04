%clear all; close all;

%% Select system to solve
file_names.starting_point = 'data/Rossler_1D_wave_train_c2.mat';
file_names.out_name = 'Rossler_spectrum_spatial_c2';

file_names.cont_func = 'absolute_spectra_rossler';
file_names.spatial_evals_fcn = 'spatial_evals_fcn_rossler';

%%
soln = load(file_names.starting_point);
u_infty = soln.U;
numPar = soln.numPar;
par = soln.par;
par.Lx = 2*pi;

%%
evals_temp = [-0.1, -0.05, 0.1];

%%
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',5000);
addpath ../utilities/
fh_spatial = str2func(file_names.spatial_evals_fcn);
[L1, L2] = ComputeLinearOperator_1D(par,numPar);

%% Find spatial eigenvalues
for k=1:length(evals_temp)
	par.lambda = evals_temp(k);
	[A,vals] = fh_spatial(L1,L2,u_infty,par,numPar);
	vals = sort(vals,'ComparisonMethod','real');
	vals = [real(vals), imag(vals)];
	save([file_names.out_name,'_',num2str(k),'.txt'],'vals','-ascii');
end

%% Plot spatial spectra
% par.lambda = 0.1;
% [A,vals] = fh_spatial(L1,L2,u_infty,par,numPar);
% vals = sort(vals,'ComparisonMethod','real');
% vals(end/2)
% vals(end/2+1)
%
% figure(1);
% clf(1);
% hold('on');
% plot(vals(1:end/2),'b.','MarkerSize',28);
% plot(vals(end/2+1:end),'r.','MarkerSize',28);
% hold off;
% xlabel('Re(\nu)'); ylabel('Im(\nu)'); box on;
% xlim([-2 2]);
% grid on;
