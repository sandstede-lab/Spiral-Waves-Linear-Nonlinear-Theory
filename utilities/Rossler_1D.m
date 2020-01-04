function [f,J]  = Rossler_1D(u,L1,L2,par,numPar,phase_cond)

nx = numPar.nx;
ny = numPar.ny;

par.(par.Free) = u(end);
a = par.a;
b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;
kappa = par.kappa;
omega = par.omega;

U = u(1:nx*ny);
V = u(nx*ny+1:2*nx*ny);
W = u(2*nx*ny+1:3*nx*ny);

% Derivative of phase condition
u_phase_th = L1*phase_cond.u_old; 

% Reaction terms
fU = -V - W;
fV = U + a*V;
fW = U.*W - c.*W + b;

line1 = delta1.*kappa^2.*(L2*U) + omega.*(L1*U) + fU;
line2 = delta2.*kappa^2.*(L2*V) + omega.*(L1*V) + fV;
line3 = delta3.*kappa^2.*(L2*W) + omega.*(L1*W) + fW;
line4 = (u_phase_th)'*(U - phase_cond.u_old);

f = [line1; line2; line3; line4];
   
% Jacobian              
if (nargout > 1)
    
fU_V = -1;  fU_W = -1;
fV_U = 1;   fV_V = a;
fW_U = W;   fW_W = U - c;
    
I = speye(nx*ny,nx*ny);
    
dU = [delta1.*kappa^2.*L2 + omega.*L1;
	 fV_U.*I;
	 spdiags(fW_U,0,nx*ny,nx*ny);
	 u_phase_th'];
    
dV = [fU_V.*I;
	  delta2.*kappa^2.*L2 + omega.*L1 + fV_V.*I;
	  sparse(nx*ny,nx*ny);
	  sparse(1,nx*ny)];
    
dW = [fU_W.*I;
	  sparse(nx*ny,nx*ny);
	  delta3.*kappa^2.*L2 + omega.*L1 + spdiags(fW_W,0,nx*ny,nx*ny);
	  sparse(1,nx*ny)];

h = 1.0e-06;
v = u;
v(end) = u(end)+h;
fp = Rossler_1D(v,L1,L2,par,numPar,phase_cond);
v(end) = u(end)-h;
fm = Rossler_1D(v,L1,L2,par,numPar,phase_cond);
dpar = (fp-fm)/(2*h);

J = [dU, dV, dW, dpar];
J = sparse(J);
             
end
