function [f,J]  = jacobian_for_spectra_fhny(u,L1,L2,par,numPar)
% J = Jacobian on the short grid format
% Need u, L1, L2 in short grid format

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
ep = par.ep;
delta = par.delta;  % Diffusion coefficient of v-equation
omega = par.omega;

%R = par.r2;
%L2 = L2./(R.^2); % Scale the laplacian by the outer radius

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);

% FHN model
fU = (1/ep)*(U.*(1-U).*(U-0.5) - (V+b)/a);
line1 =       L2*U + omega*L1*U + fU;
line2 = delta*L2*V + omega*L1*V + U - V;
f = [line1;line2];

% Jacobian              
fU_U = (1/ep)*( (1-U).*(U-0.5) - U.*(U-0.5) + U.*(1-U) );
fU_V = (1/ep)*(-1/a);
fV_U = 1;
fV_V = -1;

I = speye(nx*(ny-1)+1,nx*(ny-1)+1);
dU = [L2 + omega*L1 + spdiags(fU_U,0,nx*(ny-1)+1,nx*(ny-1)+1); fV_U.*I];
dV = [fU_V.*I; delta*L2 + omega*L1 + fV_V.*I];
J = [dU,dV];
