function [f,J]  = Karma_1D(u,L1,L2,par,numPar,phase_cond)
% 1D wave train in 2pi  frame with free parameter omega

par.(par.Free) = u(end);
tauE = 1/par.tauE;
taun = 1/par.taun;
gamma = par.gamma;
kappa = par.kappa;
delta = par.delta;
omega = par.omega;

R = 1/(1 - exp(-par.Re));

E = u(1:numPar.nx*numPar.ny);
n = u(numPar.nx*numPar.ny+1:2*numPar.nx*numPar.ny);

% Derivative of phase condition
u_phase_th = L1*phase_cond.u_old; 

% Karma model
[thSn,thSn_prime] = thetaS(E - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(E - par.Eh);        % Only calculate once

fE = tauE .* ( -E + 0.5*(par.Estar - n.^par.M).*(1 - tanE) .* E.^2 );   % Nonlinear terms
fn = taun .* ( R .* thSn - n );

line1 = kappa^2*gamma*(L2*E) + omega*(L1*E) + fE;
line2 = kappa^2*delta*(L2*n) + omega*(L1*n) + fn;
line3 = (u_phase_th)'*(E - phase_cond.u_old);       % Phase condition

f = [line1; line2; line3];
   
% Jacobian              
if (nargout > 1)
    
fE_E = tauE*(-1 - 0.5.*(par.Estar - n.^par.M).*(sech(par.Eh - E).^2).* (E.^2) + (par.Estar - n.^par.M).*(1 - tanE).*E);
fE_n = tauE*(-0.5* par.M.*n.^(par.M-1) .* (1 - tanE).*(E.^2));
fn_E = taun*(R.*thSn_prime);
fn_n = -taun;

   dE = [kappa^2*gamma*L2 + omega*L1 + spdiags(fE_E ,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny);
            spdiags(fn_E,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny) ;
            u_phase_th'];
            
   dn = [spdiags(fE_n,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny);
       kappa^2*delta*L2 + omega*L1 + fn_n*speye(numPar.nx*numPar.ny);
       sparse(1,numPar.nx*numPar.ny)];

h = 1.0e-06;
v = u;
v(end) = u(end)+h;
fp = Karma_1D(v,L1,L2,par,numPar,phase_cond);
v(end) = u(end)-h;
fm = Karma_1D(v,L1,L2,par,numPar,phase_cond);
dpar = (fp-fm)/(2*h);

J = [dE, dn, dpar];
J = sparse(J);
          
end
