function [f,J] = absolute_spectra_rossler(u,u_infty,L1,L2,par,numPar,phase_cond)
% absolute spectra of Rossler Model

nx = numPar.nx;

a = par.a;
%b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;
kappa = par.kappa;
omega = par.omega;
beta = par.beta;

U = u(1:nx); % Eigenvectors 
V = u(nx+1:2*nx);
W = u(2*nx+1:3*nx);

W1 = u(3*nx+1:4*nx);    % Second set of eigenvectors (just the imaginary parts). 
W2 = u(4*nx+1:5*nx);
W3 = u(5*nx+1:6*nx);
nu = u(end-1);          % Free parameters
lambda = u(end);

U_infty = u_infty(1:nx); % 1D asymptotic wave train
%V_infty = u_infty(nx+1:2*nx);
W_infty = u_infty(2*nx+1:3*nx);

% Phase Condtiion
u_phase = phase_cond.u_old; 
w_phase = phase_cond.w_old;

% Linearizations of non-linear terms around the wave trains
 fU_V = -1;  fU_W = -1;
 fV_U = 1;   fV_V = a;
 fW_U = W_infty;   fW_W = U_infty - c;

line1 = kappa^2*delta1*(L2*U) + 2*delta1*kappa*nu*(L1*U) + delta1*nu^2*U + omega*(L1*U) + fU_V.*V + fU_W.*W - lambda*U;    % D1(lambda,nu)(E,n) = 0
line2 = kappa^2*delta2*(L2*V) + 2*delta2*kappa*nu*(L1*V) + delta2*nu^2*V + omega*(L1*V) + fV_U.*U + fV_V.*V - lambda*V;    % D2(lambda,nu)(E,n) = 0
line3 = kappa^2*delta3*(L2*W) + 2*delta3*kappa*nu*(L1*W) + delta3*nu^2*W + omega*(L1*W) + fW_U.*U + fW_W.*W - lambda*W;    % D2(lambda,nu)(E,n) = 0


line4 = 2*delta1*kappa*(L1*U) + 2*delta1*nu*U + 1i*delta1*beta*U + 2*1i*delta1*kappa*beta*(L1*W1) + 2*1i*delta1*nu*beta*W1 - delta1*beta^2*W1 ...
    + kappa^2*delta1*(L2*W1) + 2*delta1*kappa*nu*(L1*W1) + delta1*nu^2*W1 + omega*(L1*W1) + fU_V.*W2 + fU_W.*W3 - lambda*W1; % delta1*[2(kappa*L1 + nu) + 1i*beta]*(E + 1i*beta*W1) + D1(lambda,nu)(W1,W2) = 0

line5 = 2*delta2*kappa*(L1*V) + 2*delta2*nu*V + 1i*delta2*beta*V + 2*1i*delta2*beta*kappa*(L1*W2) + 2*1i*delta2*nu*beta*W2 - delta2*beta^2*W2 ...
    + kappa^2*delta2*(L2*W2) + 2*delta2*kappa*nu*(L1*W2) + delta2*nu^2*W2 + omega*(L1*W2) + fV_U.*W1 + fV_V.*W2 - lambda*W2; % delta*[2(kappa*L1 + nu) + 1i*beta]*(n + 1i*beta*W2) + D2(lambda,nu)(W1,W2) = 0

line6 = 2*delta3*kappa*(L1*W) + 2*delta3*nu*W + 1i*delta3*beta*W + 2*1i*delta3*beta*kappa*(L1*W3) + 2*1i*delta3*nu*beta*W3 - delta3*beta^2*W3 ...
    + kappa^2*delta3*(L2*W3) + 2*delta3*kappa*nu*(L1*W3) + delta3*nu^2*W3 + omega*(L1*W3) + fW_U.*W1 + fW_W.*W3 - lambda*W3; % delta*[2(kappa*L1 + nu) + 1i*beta]*(n + 1i*beta*W2) + D2(lambda,nu)(W1,W2) = 0


line7 = (u_phase'*[U;V;W]) - 1;  % Phase condition 1
line8 = (w_phase'*[U;V;W]) - (u_phase'*[W1;W2;W3]) - 1i*beta.*(w_phase'*[W1;W2;W3]);

f = [line1; line2; line3; line4; line5; line6; line7; line8];

% Jacobian              
if (nargout > 1)
    I = speye(nx,nx); Z = sparse(nx,nx);

    
dU1 = kappa^2*delta1*L2 + 2*delta1*kappa*nu*L1 + omega*L1 +(delta1*nu^2  - lambda).*I;
dU2 = fV_U.*I; 
dU3 = spdiags(fW_U,0,nx,nx);
dU4 = 2*delta1*kappa*L1 + (2*delta1*nu + 1i*delta1*beta).*I;
dU5 = Z; 
dU6 = Z; 
dU7 = u_phase(1:nx)';
dU8 = w_phase(1:nx)';
dU = [dU1; dU2; dU3; dU4; dU5; dU6; dU7; dU8];

dV = [fU_V.*I;
    kappa^2*delta2*L2 + 2*delta2*kappa*nu*L1 + omega*L1 + (delta2*nu^2 + fV_V - lambda).*I;
    Z;Z;
    2*delta2*kappa*L1 + (2*delta2*nu + 1i*delta2*beta).*I;
    Z; u_phase(nx+1:2*nx)';w_phase(nx+1:2*nx)'];

dW = [fU_W.*I;Z;
    kappa^2*delta3*L2 + 2*delta3*kappa*nu*L1 + omega*L1 + spdiags(delta3*nu^2 + fW_W - lambda,0,nx,nx);
    Z;Z;
    2*delta3*kappa*L1 + (2*delta3*nu + 1i*delta3*beta).*I;
    u_phase(2*nx+1:3*nx)';w_phase(2*nx+1:3*nx)'];

dW1 = [Z;Z;Z;
    2*1i*delta1*kappa*beta*L1 + kappa^2*delta1*L2 + 2*delta1*kappa*nu*L1 + omega*L1 + (2*1i*delta1*nu*beta - delta1*beta^2 + delta1*nu^2  - lambda).*I;
    fV_U.*I; spdiags(fW_U,0,nx,nx);
    sparse(1,nx);
    -u_phase(1:nx)' - 1i*beta*w_phase(1:nx)'];
    

dW2 = [ Z; Z; Z;
    fU_V.*I;
    2*1i*delta2*beta*kappa*L1 + kappa^2*delta2*L2 + 2*delta2*kappa*nu*L1 + omega*L1 + (2*1i*delta2*nu*beta - delta2*beta^2 + delta2*nu^2 + fV_V - lambda)*speye(nx);
    Z; sparse(1,nx);
    -u_phase(nx+1:2*nx)' - 1i*beta*w_phase(nx+1:2*nx)'];

dW3 = [ Z; Z; Z;
    fU_W.*I;
    Z;
    2*1i*delta3*beta*kappa*L1 + kappa^2*delta3*L2 + 2*delta3*kappa*nu*L1 + omega*L1 + spdiags(2*1i*delta3*nu*beta - delta3*beta^2 + delta3*nu^2 + fW_W - lambda,0,nx,nx);
    sparse(1,nx);
    -u_phase(2*nx+1:3*nx)' - 1i*beta*w_phase(2*nx+1:3*nx)'];

dnu = [ 2*delta1*kappa*(L1*U) + 2*delta1*nu*U;
    2*delta2*kappa*(L1*V) + 2*delta2*nu*V;
    2*delta3*kappa*(L1*W) + 2*delta3*nu*W;
    2*delta1*U + 2*1i*delta1*beta*W1 + 2*delta1*kappa*(L1*W1) + 2*delta1*nu*W1;
    2*delta2*V + 2*1i*delta2*beta*W2 + 2*delta2*kappa*(L1*W2) + 2*delta2*nu*W2;
    2*delta3*W + 2*1i*delta3*beta*W3 + 2*delta3*kappa*(L1*W3) + 2*delta3*nu*W3;
    0;
    0];

dlambda = [-U;-V;-W; -W1; -W2; -W3; 0; 0];

J = [dU, dV, dW, dW1, dW2, dW3, dnu, dlambda];
   
J = sparse(J);
end







