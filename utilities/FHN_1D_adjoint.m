function J  = FHN_1D_adjoint(u,L1,L2,par,numPar)
% Asymptotics of FHN model with diffusion in the v-equation

nx = numPar.nx;

a = par.a;
b = par.b;
ep = par.ep;
kappa = par.kappa;
delta = par.delta; % v-diffusion coefficient
omega = par.omega;

U = u(1:nx);
V = u(nx+1:2*nx);

fU_U = (1/ep)*( (1-U).*(U-0.5) - U.*(U-0.5) + U.*(1-U) );
fU_V = (1/ep)*(-1/a);
I = speye(nx,nx);

dU = [kappa^2*L2 + omega*L1 + spdiags(fU_U,0,nx,nx); I];
dV = [fU_V*I; delta*kappa^2*L2 + omega*L1 - I];

J = [dU, dV];
J = sparse(J);

end
