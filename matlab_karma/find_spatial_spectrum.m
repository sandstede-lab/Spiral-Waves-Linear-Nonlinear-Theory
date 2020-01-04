%clear all; close all;

%% Select system to solve
file_names.starting_point = 'Karma_1D_wave_train_re1p2.mat';
file_names.out_name = 'Karma_spatial_spectrum.mat';
file_names.spatial_evals_fcn = 'spatial_evals_fcn_karma';

%%
soln = load(file_names.starting_point);
u_infty = soln.U;
numPar = soln.numPar;
par = soln.par;
par.Lx = 2*pi;

%%
%par.lambda = -12 + i*par.omega/2 + 1.7*i;
%par.lambda = -2.20 + 1*i*28.12 + 3.7;
par.lambda = -9.4;

%%
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
	'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',5000);

addpath ../utilities/
fh_spatial = str2func(file_names.spatial_evals_fcn);

%%
[L1, L2] = ComputeLinearOperator_1D(par,numPar);

%% Find spatial eigenvalues
[A,V,vals] = fh_spatial(L1,L2,u_infty,par,numPar);
[vals, idx] = sort(vals,'ComparisonMethod','real');
V = V(:,idx);

vals(end/2)
vals(end/2+1)

%% Plot spatial eigenvalues
figure(1);
clf(1);
hold('on');
plot(vals(1:end/2),'b.','MarkerSize',28);
plot(vals(end/2+1:end),'r.','MarkerSize',28);
hold off;
xlabel('Re(\nu)'); ylabel('Im(\nu)'); box on;
xlim([-20 20]);
grid on;

%% absolute spectrum
%par.nu = vals(end/2);
%par.beta = imag(vals(end/2+1)-vals(end/2));
%uout = V(:,end/2);

%% essential spectrum
nx = numPar.nx;
k = 0;
par.nu = vals(end/2+k);
uout = [V(1:nx,end/2+k); V(2*nx+1:3*nx,end/2+k)];

%% save data
save(file_names.out_name ,'uout','par','numPar','u_infty');
