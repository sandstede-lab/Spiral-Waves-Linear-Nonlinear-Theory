function [f,J]  = Karma_essSpec_fcn(u,u_infty,L1,L2,par,numPar,phase_cond)
% Spiral dispersion relation with free parameter the temporal eigenvalue.
% Used in computation of essential spectrum
% In spiral coordinates

nx = numPar.nx;
ny = numPar.ny;

tauE = 1/par.tauE;
taun = 1/par.taun;
gamma = par.gamma;
kappa = par.kappa;
delta = par.delta;
nu = par.nu;
omega = par.omega;

R = 1/(1 - exp(-par.Re));

E = u(1:nx); % Eigenvectors 
n = u(nx+1:2*nx);
lambda = u(end); % Free parameter = temporal eigenvalue, lambda


E_infty = u_infty(1:nx); % 1D asymptotic wave train
n_infty = u_infty(nx+1:2*nx);

% Phase Condtiion
u_phase = phase_cond.u_old;

% Karma model
[~,thSn_prime] = thetaS(E_infty - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(E_infty - par.Eh); 

% Linearizations of non-linear terms around the wave trains
fE_E = tauE*(-1 - 0.5.*(par.Estar - n_infty.^par.M).*(sech(par.Eh - E_infty).^2).* (E_infty.^2) + (par.Estar - n_infty.^par.M).*(1 - tanE).*E_infty);
fE_n = tauE*(-0.5* par.M.*n_infty.^(par.M-1) .* (1 - tanE).*(E_infty.^2));
fn_E = taun*(R.*thSn_prime);
fn_n = -taun;

% spiral dispersion relation: nu is a COMPLEX number
line1 = kappa^2*gamma*(L2*E) + 2*gamma*kappa*nu*(L1*E) + gamma*nu^2*E + omega*(L1*E) + fE_E.*E + fE_n.*n - lambda*E;
line2 = kappa^2*delta*(L2*n) + 2*delta*kappa*nu*(L1*n) + delta*nu^2*n + omega*(L1*n) + fn_E.*E + fn_n.*n - lambda*n;
line3 = (u_phase'*[E;n]) - 1; 

f = [line1; line2; line3];
 
% Jacobian              
if (nargout > 1)

dE1 = kappa^2*gamma*L2 + 2*gamma*kappa*nu*L1 + omega*L1 + spdiags(gamma*nu^2 -lambda +fE_E ,0,nx*ny,nx*ny);
dE2 = spdiags( fn_E ,0,nx*ny,nx*ny) ;
dE3 = u_phase(1:nx)';
            
dE = [dE1; dE2; dE3];

dn = [spdiags(fE_n,0,nx*ny,nx*ny);
       kappa^2*delta*L2 + 2*delta*kappa*nu*L1 + omega*L1 + (delta*nu^2 - lambda + fn_n)*speye(nx*ny);
       u_phase(nx+1:2*nx)'];
   
   dlambda = [-E;
        -n;
        0];
   
   J = [dE, dn, dlambda];
   
   J = sparse(J);
end



    
    
    
    
    
    
    