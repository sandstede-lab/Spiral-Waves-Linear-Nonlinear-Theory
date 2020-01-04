function [A,V,vals] = spatial_evals_fcn_bar_eiswirth(L1,L2,u_infty,par,numPar)
% Finds the spatial eigenvalues. Equations have diffusion in v-equation

nx = numPar.nx;
lambda = par.lambda;
kappa = par.kappa;
ep = par.ep;
a = par.a;
b = par.b;
omega = par.omega;
delta = par.delta;

U_infty = u_infty(1:numPar.nx);
V_infty = u_infty(numPar.nx+1:2*numPar.nx);

fU_U = (1/ep)*(2*(1+b/a)*U_infty-b/a-(1/a)*V_infty-3*U_infty.^2+(2/a)*U_infty.*V_infty);
fU_V = (1/ep)*(-(1/a)*U_infty + (1/a)*U_infty.^2);

fV_U = -6.75.*((U_infty-1).^2 + 2.*U_infty.*(U_infty-1)).*(U_infty<=1).*(U_infty>=1/3);
fV_V = -1;

I = speye(nx,nx);

Ex = [sparse(nx,nx), I, sparse(nx,nx), sparse(nx,nx)];
Fx = -1*[ kappa^2*L2 + omega*L1 + spdiags(fU_U - lambda,0,nx,nx), 2*kappa*L1, spdiags(fU_V,0,nx,nx), sparse(nx,nx) ];
Nx = [ sparse(nx,nx), sparse(nx,nx), sparse(nx,nx), I ];
mx = -1/delta*[ spdiags(fV_U,0,nx,nx), sparse(nx,nx), delta*kappa^2*L2 + omega*L1 + (fV_V - lambda)*speye(nx,nx), 2*delta*kappa*L1];
A = [Ex; Fx; Nx; mx];

if nargout > 1
	[V,vals] = eig(full(A));
	vals = diag(vals);
end
