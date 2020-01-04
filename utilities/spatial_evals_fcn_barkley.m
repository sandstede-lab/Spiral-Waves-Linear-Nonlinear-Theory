function [A,vals] = spatial_evals_fcn_barkley(L1,L2,u_infty,par,numPar)
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


fU_U = (1/ep)*(-3*U_infty.^2 + 2.*U_infty.*(V_infty+b)./a + 2.*U_infty - (V_infty+b)./a); % f_u linearization
fU_V = 1/(ep*a).*U_infty.*(U_infty-1); % f_v lineariztion

I = speye(nx,nx);

Ex = [sparse(nx,nx), I, sparse(nx,nx), sparse(nx,nx)];
Fx = -1*[kappa^2*L2 + omega*L1 + spdiags(fU_U - lambda,0,nx,nx), 2*kappa*L1, spdiags(fU_V,0,nx,nx), sparse(nx,nx)];
Nx = [ sparse(nx,nx), sparse(nx,nx), sparse(nx,nx), I];
mx = -(1/delta)*[ I, sparse(nx,nx), delta*kappa^2*L2 + omega*L1 - (1 + lambda)*speye(nx,nx), 2*delta*kappa*L1];
A = [Ex; Fx; Nx; mx];

if nargout > 1

vals = eig(full(A));  
%vals = diag(vals); % Don't need if don't ask for eigenvectors.
end


