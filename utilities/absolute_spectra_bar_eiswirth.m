function [f,J] = absolute_spectra_barkley(u,u_infty,L1,L2,par,numPar,phase_cond)
% Continue absolute spectra of Barkley model 
% Follows continuation methods described in Rademacher et al. 2007

nx = numPar.nx;

a = par.a;
b = par.b;
ep = par.ep;
omega = par.omega;
kappa = par.kappa;
beta = par.beta;
delta = par.delta;

U = u(1:nx);            % First set of eigenvectors
V = u(nx+1:2*nx);
W1 = u(2*nx+1:3*nx);    % Second set of eigenvectors (just the imaginary parts). 
W2 = u(3*nx+1:4*nx);
nu = u(end-1);          % Free parameters
lambda = u(end);

U_infty = u_infty(1:nx); % 1D asymptotic wave train
V_infty = u_infty(nx+1:2*nx);

% Phase Condtiion
u_phase = phase_cond.u_old; 
w_phase = phase_cond.w_old;

% Linearizations
fU_U = (1/ep)*(-3*U_infty.^2 + 2.*U_infty.*(V_infty+b)./a + 2.*U_infty - (V_infty+b)./a); % f_u linearization
fU_V = 1/(ep*a).*U_infty.*(U_infty-1); % f_v lineariztion

fV_U = -6.75.*((U_infty-1).^2 + 2.*U_infty.*(U_infty-1)).*(U_infty<=1).*(U_infty>=1/3);
fV_V = -1;

line1 = nu^2*U + 2*kappa*nu*L1*U + kappa^2*L2*U + omega*L1*U + fU_U.*U + fU_V.*V - lambda.*U;
line2 = delta*(nu^2*V + 2*kappa*nu*L1*V + kappa^2*L2*V) + omega*L1*V + fV_U.*U + fV_V.*V  - lambda.*V;

line3 = 2*kappa*L1*U + 2*nu*U + 1i*beta*U + 2*1i*kappa*beta*L1*W1 + 2*1i*beta*nu*W1 - (beta^2)*W1 + omega*L1*W1...
    + nu^2*W1 + 2*kappa*nu*L1*W1 + kappa^2*L2*W1 + fU_U.*W1 + fU_V.*W2 - lambda.*W1;

line4 = 2*delta*kappa*L1*V + 2*delta*nu*V + 1i*delta*beta*V + 2*1i*delta*kappa*beta*L1*W2 + 2*1i*delta*beta*nu*W2 - delta*(beta^2)*W2 ...
    + delta*nu^2*W2 + 2*delta*kappa*nu*L1*W2 + kappa^2*delta*L2*W2 + fV_U.*W1 + fV_V.*W2 + omega*L1*W2 - lambda.*W2;

line5 = (u_phase'*[U;V]) - 1;   % Phase condition 1
line6 = (w_phase'*[U;V]) - (u_phase'*[W1;W2]) - 1i*beta.*(w_phase'*[W1;W2]);

f = [line1;
    line2;
    line3;
    line4;
    line5;
    line6];

% Jacobian              
if (nargout > 1)

dU = [omega*L1  + 2*kappa*nu*L1 + kappa^2*L2 + spdiags(fU_U + nu^2 - lambda,0,nx,nx);    
    spdiags(fV_U,0,nx,nx);
    2*kappa*L1 + (2*nu + 1i*beta)*speye(nx);
    zeros(nx);
    u_phase(1:nx)';
    w_phase(1:nx)'];

dV = [spdiags(fU_V,0,nx,nx);
     2*delta*kappa*nu*L1 + kappa^2*delta*L2 + omega*L1 + (delta*nu^2 - lambda)*speye(nx) + fV_V*speye(nx);
    zeros(nx);
     2*delta*kappa*L1 + (2*delta*nu + 1i*delta*beta) *speye(nx);
    u_phase(nx+1:2*nx)';
    w_phase(nx+1:2*nx)'];

dW1 = [zeros(nx);
    zeros(nx);
    2*1i*kappa*beta*L1 + omega*L1  + 2*kappa*nu*L1 + kappa^2*L2 + spdiags(fU_U + 2*1i*beta*nu - (beta^2) + nu^2 - lambda, 0, nx,nx);
    spdiags(fV_U,0,nx,nx);
    zeros(1,nx);
    -1*u_phase(1:nx)' - 1i*beta.*w_phase(1:nx)'];

dW2 = [zeros(nx);
    zeros(nx);
    spdiags(fU_V,0,nx,nx);
    2*1i*delta*kappa*beta*L1 + 2*delta*kappa*nu*L1 + kappa^2*delta*L2 + omega*L1 + fV_V*speye(nx) + (delta*nu^2- lambda+ 2*1i*delta*beta*nu - delta*(beta^2))*speye(nx);
    zeros(1,nx);
     -1*u_phase(nx+1:2*nx)' - 1i*beta.*w_phase(nx+1:2*nx)'];

 dnu = [2*nu*U + 2*kappa*L1*U;     
    2*delta*nu*V + 2*delta*kappa*L1*V;
    2*U + 2*1i*beta*W1 + 2*nu*W1 + 2*kappa*L1*W1;
    2*delta*V + 2*1i*delta*beta*W2 + 2*delta*nu*W2 + 2*delta*kappa*L1*W2;
    0;
    0];


dlambda = [-U;    
    -V;
    -W1;
    -W2;
    0;
    0]; 
J = [dU, dV, dW1, dW2, dnu, dlambda];

J = sparse(J);
end








