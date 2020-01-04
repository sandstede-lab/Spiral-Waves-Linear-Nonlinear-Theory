function [f,J]  = jacobian_for_spectra_bar_eiswirth(u,L1,L2,par,numPar)
% J = Jacobian on the short grid format
% Need u, L1, L2 in short grid format

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
ep = par.ep;
delta = par.delta;  % Diffusion coefficient of v-equation
omega = par.omega;

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);

% Bar-Eiswirth model
nl1 = (1/ep)*U.*(1-U).*(U - (V+b)/a);
nl2 = (1-6.75.*U.*(U-1).^2).*(U<=1).*(U>=1/3) + 1.*(U>1);
line1 = L2*U       + omega*L1*U + nl1;
line2 = delta*L2*V + omega*L1*V + nl2 - V;
f = [line1;line2];
   
% Jacobian              
fU_U = (1/ep)*(2*(1+b/a)*U-b/a-(1/a)*V-3*U.^2+(2/a)*U.*V);
fU_V = (1/ep)*(-(1/a)*U + (1/a)*U.^2);
fV_U = -6.75.*((U-1).^2 +2.*U.*(U-1)).*(U<=1).*(U>=1/3);
fV_V = -1;

I = speye(nx*(ny-1)+1,nx*(ny-1)+1);

dU = [L2 + omega*L1 + spdiags(fU_U,0,nx*(ny-1)+1,nx*(ny-1)+1);
	fV_U.*I];

dV = [spdiags(fU_V,0,nx*(ny-1)+1,nx*(ny-1)+1);
	delta*L2 + omega*L1 + fV_V.*I];

J = [dU,dV];
