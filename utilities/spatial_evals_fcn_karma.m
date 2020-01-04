function [A,V,vals] = spatial_evals_fcn_karma(L1,L2,u_infty,par,numPar)
% Finds the spatial eigenvalues. Equations have diffusion in v-equation

nx = numPar.nx;

E_infty = u_infty(1:numPar.nx);
n_infty = u_infty(numPar.nx+1:2*numPar.nx);

tauE = 1/par.tauE;
taun = 1/par.taun;
gamma = par.gamma;
kappa = par.kappa;
delta = par.delta;
lambda = par.lambda;
omega = par.omega;

R = 1/(1 - exp(-par.Re));

[~,thSn_prime] = thetaS(E_infty - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(E_infty - par.Eh );

% Linearizations of non-linear terms around the wave trains
fE_E = tauE*(-1 - 0.5.*(par.Estar - n_infty.^par.M).*(sech(par.Eh - E_infty).^2).* (E_infty.^2) + (par.Estar - n_infty.^par.M).*(1 - tanE).*E_infty);
fE_n = tauE*(-0.5* par.M.*n_infty.^(par.M-1) .* (1 - tanE).*(E_infty.^2));
fn_E = taun*(R.*thSn_prime);
fn_n = -taun;

Ex = [sparse(nx,nx), speye(nx,nx), sparse(nx,nx), sparse(nx,nx)];
Fx = -(1/gamma)*[gamma*kappa^2*L2 + omega*L1 + spdiags(fE_E - lambda,0,nx,nx), 2*kappa*gamma*L1, spdiags(fE_n,0,nx,nx), sparse(nx,nx)];
Nx = [ sparse(nx,nx), sparse(nx,nx), sparse(nx,nx), speye(nx,nx)];
mx = -(1/delta)*[ spdiags(fn_E,0,nx,nx), sparse(nx,nx), delta*kappa^2*L2 + omega*L1 + (fn_n - lambda)*speye(nx,nx), 2*delta*kappa*L1];

A = [Ex; Fx; Nx; mx];

if nargout > 1
	[V,vals] = eig(full(A));
	vals = diag(vals);
end


