function [A,V,vals] = spatial_evals_fcn_rossler(L1,L2,u_infty,par,numPar)
% Finds the spatial eigenvalues. Equations have diffusion in v-equation

nx = numPar.nx;

U_infty = u_infty(1:numPar.nx);
%V_infty = u_infty(numPar.nx+1:2*numPar.nx);
W_infty = u_infty(2*numPar.nx+1:3*numPar.nx);

a = par.a;
%b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;
kappa = par.kappa;
omega = par.omega;
lambda = par.lambda;


% Linearizations of non-linear terms around the wave trains
 fU_V = -1;  fU_W = -1;
 fV_U = 1;   fV_V = a;
 fW_U = W_infty;   fW_W = U_infty - c;

 % Need to work on this part! - write out how the 3 matrix works!
 Z = sparse(nx,nx); I = speye(nx,nx);
 
Ux = [ Z, I, Z, Z, Z, Z];
Ubarx = -(1/delta1).*[delta1.*kappa^2.*L2 + (omega*L1 + - lambda.*I), 2*kappa*delta1*L1, fU_V.*I, Z, fU_W.*I, Z];
Vx = [Z, Z, Z, I, Z, Z];
Vbarx = -(1/delta2).*[fV_U.*I, Z, kappa^2.*delta2.*L2 + (omega*L1 +(fV_V - lambda).*I), 2*kappa*delta2*L1, Z, Z ];
Wx = [Z, Z, Z, Z, Z, I];
Wbarx = -(1/delta3).*[ spdiags(fW_U,0,nx,nx), Z, Z, Z, delta3.*kappa^2.*L2 + (omega*L1 + spdiags(fW_W - lambda,0,nx,nx)), 2*kappa*delta3*L1  ];

A = [Ux; Ubarx; Vx; Vbarx; Wx; Wbarx];

if nargout > 1
%numEvals = 8;
%numEvals = 256;
%vals = eigs(A,numEvals,0);  % Do I really trust this?
[V,vals] = eig(full(A));

vals = diag(vals); % Don't need if don't ask for eigenvectors.
end
