function [f,J]  = FHN_1D(u,L1,L2,par,numPar,phase_cond)
% Asymptotics of FHN model with diffusion in the v-equation

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

%% phase condition
u_phase_th = L1*phase_cond.u_old; 

%% FHN model
fU = (1/ep)*(U.*(1-U).*(U-0.5) - (V+b)/a);
f = [      kappa^2*(L2*U) + omega*L1*U + fU;
	 delta*kappa^2*(L2*V) + omega*L1*V + U-V ;  
	 (u_phase_th)'*(U - phase_cond.u_old)];

%% Jacobian              
if (nargout > 1)
	fU_U = (1/ep)*( (1-U).*(U-0.5) - U.*(U-0.5) + U.*(1-U) );
	fU_V = (1/ep)*(-1/a);
	I = speye(nx,nx);

	dU = [kappa^2*L2 + omega*L1 + spdiags(fU_U,0,nx,nx); I; u_phase_th'];
	dV = [fU_V*I; delta*kappa^2*L2 + omega*L1 - I; sparse(1,nx)];

	h = 1.0e-06;
	v = u;
	v(end) = u(end)+h;
	fp = FHN_1D(v,L1,L2,par,numPar,phase_cond);
	v(end) = u(end)-h;
	fm = FHN_1D(v,L1,L2,par,numPar,phase_cond);
	dpar = (fp-fm)/(2*h);

	J = [dU, dV, dpar];
	J = sparse(J);    
end
