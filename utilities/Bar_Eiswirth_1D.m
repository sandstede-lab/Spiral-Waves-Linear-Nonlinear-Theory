function [f,J]  = Bar_Eiswirth_1D(u,L1,L2,par,numPar,phase_cond)
% Asymptotics of Bar_Eiswirth model with diffusion in the v-equation

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

% Bar-Eiswirth model
nl1 = (1/ep)*U.*(1-U).*(U - (V+b)/a);
nl2 = (1-6.75.*U.*(U-1).^2).*(U<=1).*(U>=1/3) + 1.*(U>1);
f = [omega*(L1*U) + kappa^2*(L2*U) + nl1;
	 delta*kappa^2*(L2*V) + nl2 - V + omega*L1*V;  
	 (u_phase_th)'*(U - phase_cond.u_old)];

% Jacobian              
if (nargout > 1)
    fU_U = (1/ep)*(2*(1+b/a)*U-b/a-(1/a)*V-3*U.^2+(2/a)*U.*V);
    fU_V = (1/ep)*(-(1/a)*U + (1/a)*U.^2);

    fV_U = -6.75.*((U-1).^2 + 2.*U.*(U-1)).*(U<=1).*(U>=1/3);
    fV_V = -1;
	I = speye(nx,nx);

	% dU, dV, dpar

	dU = [kappa^2*L2 + omega*L1 + spdiags(fU_U,0,nx,nx);
		fV_U.*I; u_phase_th'];

	dV1 = spdiags(fU_V,0,nx,nx);
	dV2 = delta*kappa^2*L2 +omega*L1 + fV_V.*I;
	dV3 = sparse(1,nx);
	dV = [dV1; dV2; dV3];

	h = 1.0e-06;
	v = u;
	v(end) = u(end)+h;
	fp = Bar_Eiswirth_1D(v,L1,L2,par,numPar,phase_cond);
	v(end) = u(end)-h;
	fm = Bar_Eiswirth_1D(v,L1,L2,par,numPar,phase_cond);
	dpar = (fp-fm)/(2*h);

	J = [dU, dV, dpar];
	J = sparse(J);
end
