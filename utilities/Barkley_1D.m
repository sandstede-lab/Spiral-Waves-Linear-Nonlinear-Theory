function [f,J]  = Barkley_1D(u,L1,L2,par,numPar,phase_cond)
% Asymptotics of Barkley model with diffusion in the v-equation

nx = numPar.nx;

par.(par.Free) = u(end);
a = par.a;
b = par.b;
ep = par.ep;
kappa = par.kappa;
delta = par.delta; % v-diffusion coefficient
omega = par.omega;

U = u(1:nx);
V = u(nx+1:2*nx);

% phase condition
u_phase_th = L1*phase_cond.u_old; 

% Barkley model
nonlin = (1/ep)*U.*(1-U).*(U - (V+b)/a);
f = [omega*(L1*U) + kappa^2*(L2*U) + nonlin;
	 delta*kappa^2*(L2*V) + U-V + omega*L1*V;  
	 (u_phase_th)'*(U - phase_cond.u_old)];

% Jacobian              
if (nargout > 1)
nonlin_jacob1 = (1/ep)*(2*(1+b/a)*U-b/a-(1/a)*V-3*U.^2+(2/a)*U.*V);
nonlin_jacob2 = (1/ep)*(-(1/a)*U + (1/a)*U.^2);
I = speye(nx,nx);

% dU, dV, dpar

dU = [kappa^2*L2 + omega*L1 + spdiags(nonlin_jacob1,0,nx,nx); I; u_phase_th'];

dV1 = spdiags(nonlin_jacob2,0,nx,nx);
dV2 = delta*kappa^2*L2 +omega*L1 - I;
dV3 = sparse(1,nx);
dV = [dV1; dV2; dV3];

h = 1.0e-06;
v = u;
v(end) = u(end)+h;
fp = Barkley_1D(v,L1,L2,par,numPar,phase_cond);
v(end) = u(end)-h;
fm = Barkley_1D(v,L1,L2,par,numPar,phase_cond);
dpar = (fp-fm)/(2*h);

J = [dU, dV, dpar];
J = sparse(J);
          
end
