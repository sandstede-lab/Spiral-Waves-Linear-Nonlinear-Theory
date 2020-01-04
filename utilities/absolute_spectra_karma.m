function [f,J] = absolute_spectra_karma(u,u_infty,L1,L2,par,numPar,phase_cond)
% Continue absolute spectra of Karma model (fromKarma 1994)
% Follows continuation methods described in Rademacher et al. 2007

nx = numPar.nx;
ny = numPar.ny;

tauE = 1/par.tauE;
taun = 1/par.taun;
gamma = par.gamma;
kappa = par.kappa;
delta = par.delta;
omega = par.omega;
En = par.En;
Eh = par.Eh;
s = par.s;

beta = par.beta; % Continuation parameter

R = 1 - exp(-par.Re);

E = u(1:nx); % Eigenvectors 
n = u(nx+1:2*nx);
W1 = u(2*nx+1:3*nx);    % Second set of eigenvectors (just the imaginary parts). 
W2 = u(3*nx+1:4*nx);
nu = u(end-1);          % Free parameters
lambda = u(end);

E_infty = u_infty(1:nx); % 1D asymptotic wave train
n_infty = u_infty(nx+1:2*nx);

% Phase Condtiion
u_phase = phase_cond.u_old; 
w_phase = phase_cond.w_old;

% Karma model
[thSn,thSn_prime] = thetaS(E_infty - En,s);   % Smoothed Heaviside function
tanE = tanh(Eh - E_infty );

% Linearizations of non-linear terms around the wave trains
fE_E = tauE*(-1 - 0.5.*(par.Estar - n_infty.^par.M).*(sech(par.Eh - E_infty).^2).* (E_infty.^2) + (par.Estar - n_infty.^par.M).*(1 + tanE).*E_infty);
fE_n = tauE*(-0.5* par.M.*n_infty.^(par.M-1) .* (1 + tanE).*(E_infty.^2));
fn_E = taun*(((1 - R.* n_infty)./R).*thSn_prime  + n_infty.*thSn_prime);
fn_n = -taun;

line1 = kappa^2*gamma*(L2*E) + 2*gamma*kappa*nu*(L1*E) + gamma*nu^2*E + omega*(L1*E) + fE_E.*E + fE_n.*n - lambda*E;    % D1(lambda,nu)(E,n) = 0
line2 = kappa^2*delta*(L2*n) + 2*delta*kappa*nu*(L1*n) + delta*nu^2*n + omega*(L1*n) + fn_E.*E + fn_n.*n - lambda*n;    % D2(lambda,nu)(E,n) = 0


line3 = 2*gamma*kappa*(L1*E) + 2*gamma*nu*E + 1i*gamma*beta*E + 2*1i*gamma*kappa*beta*(L1*W1) + 2*1i*gamma*nu*beta*W1 - gamma*beta^2*W1 ...
    + kappa^2*gamma*(L2*W1) + 2*gamma*kappa*nu*(L1*W1) + gamma*nu^2*W1 + omega*(L1*W1) + fE_E.*W1 + fE_n.*W2 - lambda*W1; % gamma*[2(kappa*L1 + nu) + 1i*beta]*(E + 1i*beta*W1) + D1(lambda,nu)(W1,W2) = 0

line4 = 2*delta*kappa*(L1*n) + 2*delta*nu*n + 1i*delta*beta*n + 2*1i*delta*beta*kappa*(L1*W2) + 2*1i*delta*nu*beta*W2 - delta*beta^2*W2 ...
    + kappa^2*delta*(L2*W2) + 2*delta*kappa*nu*(L1*W2) + delta*nu^2*W2 + omega*(L1*W2) + fn_E.*W1 + fn_n.*W2 - lambda*W2; % delta*[2(kappa*L1 + nu) + 1i*beta]*(n + 1i*beta*W2) + D2(lambda,nu)(W1,W2) = 0

line5 = (u_phase'*[E;n]) - 1;  % Phase condition 1
line6 = (w_phase'*[E;n]) - (u_phase'*[W1;W2]) - 1i*beta.*(w_phase'*[W1;W2]);

f = [line1;
    line2;
    line3;
    line4;
    line5;
    line6];

% Jacobian              
if (nargout > 1)

dE = [kappa^2*gamma*L2 + 2*gamma*kappa*nu*L1 + omega*L1 + spdiags(gamma*nu^2 + fE_E - lambda, 0, nx,nx);
    spdiags(fn_E,0,nx,nx);
    2*gamma*kappa*L1 + (2*gamma*nu + 1i*gamma*beta)*speye(nx);
    sparse(nx,nx);
    u_phase(1:nx)';
    w_phase(1:nx)'];

dn = [spdiags(fE_n,0,nx,nx);
    kappa^2*delta*L2 + 2*delta*kappa*nu*L1 + omega*L1 + (delta*nu^2 + fn_n - lambda)*speye(nx);
    sparse(nx,nx);
    2*delta*kappa*L1 + (2*delta*nu + 1i*delta*beta)*speye(nx);
    u_phase(nx+1:2*nx)';
    w_phase(nx+1:2*nx)'];

dW1 = [sparse(nx,nx);
    sparse(nx,nx);
    2*1i*gamma*kappa*beta*L1 + kappa^2*gamma*L2 + 2*gamma*kappa*nu*L1 + omega*L1 + spdiags(2*1i*gamma*nu*beta - gamma*beta^2 + gamma*nu^2 + fE_E - lambda,0,nx,nx);
    spdiags(fn_E,0,nx,nx);
    sparse(1,nx);
    -u_phase(1:nx)' - 1i*beta*w_phase(1:nx)'];
    

dW2 = [ sparse(nx,nx);
    sparse(nx,nx);
    spdiags(fE_n,0,nx,nx);
    2*1i*delta*beta*kappa*L1 + kappa^2*delta*L2 + 2*delta*kappa*nu*L1 + omega*L1 + (2*1i*delta*nu*beta - delta*beta^2 + delta*nu^2 + fn_n - lambda)*speye(nx);
    sparse(1,nx);
    -u_phase(nx+1:2*nx)' - 1i*beta*w_phase(nx+1:2*nx)'];

dnu = [ 2*gamma*kappa*(L1*E) + 2*gamma*nu*E;
    2*delta*kappa*(L1*n) + 2*delta*nu*n;
    2*gamma*E + 2*1i*gamma*beta*W1 + 2*gamma*kappa*(L1*W1) + 2*gamma*nu*W1;
    2*delta*n + 2*1i*delta*beta*W2 + 2*delta*kappa*(L1*W2) + 2*delta*nu*W2;
    0;
    0];

dlambda = [-E;
    -n;
    -W1;
    -W2;
    0;
    0];

J = [dE, dn, dW1, dW2, dnu, dlambda];
   
J = sparse(J);
end







