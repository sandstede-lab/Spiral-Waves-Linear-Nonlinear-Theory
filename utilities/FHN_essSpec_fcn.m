function [f,J] = FHN_essSpec_fcn(u,u_infty,L1,L2,par,numPar,phase_cond)
% Spiral dispersion relation with free parameter the temporal eigenvalue.
% Used in computation of essential spectrum
% In spiral coordinates

nx = numPar.nx;

a = par.a;
b = par.b;
ep = par.ep;
omega = par.omega;
nu = par.nu;
kappa = par.kappa;
delta = par.delta;
flag = par.flag;	% =0 for transverse; =1 for normal essential spectrum

U = u(1:nx); % Eigenvectors 
V = u(nx+1:2*nx);
lambda = u(end); % Free parameter = temporal eigenvalue, lambda

U_infty = u_infty(1:nx); % 1D asymptotic wave train
V_infty = u_infty(nx+1:2*nx);

%% Phase Condition
u_phase = phase_cond.u_old; 

%% FHN model
fU_U = (1/ep)*( (U_infty-0.5).*(1-U_infty) + U_infty.*(1-U_infty) - U_infty.*(U_infty-0.5) );
fU_V = (1/ep)*(-1/a);

line1 =        kappa^2*L2*U + flag*2*kappa*nu*L1*U + nu^2*U  + omega*L1*U + fU_U.*U + fU_V*V - lambda*U;
line2 = delta*(kappa^2*L2*V + flag*2*kappa*nu*L1*V + nu^2*V) + omega*L1*V + U - V - lambda*V;
line3 = (u_phase'*[U;V]) - 1; % Phase condition
f = [line1; line2; line3];
   
%% Jacobian              
if (nargout > 1)
	I = speye(nx,nx);
	dU = [ kappa^2*L2 + flag*2*kappa*nu*L1 + nu^2*I + omega*L1 + spdiags(fU_U-lambda,0,nx,nx);
		I;
		u_phase(1:nx)'];
	dV = [fU_V*I;
		delta*(kappa^2*L2 + flag*2*kappa*nu*L1 + nu^2*I) + omega*L1 - (1+lambda)*I;
		u_phase(nx+1:2*nx)'];
	dlambda = [-U; -V; 0];

	J = [dU, dV, dlambda];  
	J = sparse(J);
end
