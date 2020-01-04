function [f,J] = jacobian_for_spectra_karma(u,L1,L2,par,numPar)
% J = Jacobian on the short grid format
% Need u, L1, L2 in short grid format

% Parameters
tauE = 1/par.tauE;
taun = 1/par.taun;

Re = 1/(1 - exp(-par.Re));
omega = par.omega;

nx = numPar.nx;
ny = numPar.ny;

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);

[thSn,thSn_prime] = thetaS(U - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(U - par.Eh);        % Only calculate once

% Karma model

fU = tauE .* ( -U + 0.5*(par.Estar - V.^par.M).*(1 - tanE) .* U.^2 );
fV = taun .* ( Re.* thSn - V );

line1 = par.gamma*(L2*U) + omega*(L1*U) + fU;
line2 = par.delta*(L2*V) + omega*(L1*V) + fV;

f = [line1; line2];

fU_U = tauE*(-1 - 0.5.*(par.Estar - V.^par.M).*(sech(par.Eh - U).^2).* (U.^2) + (par.Estar - V.^par.M).*(1 - tanE).*U);  % This is correct
fU_V = tauE*(-0.5* par.M.*V.^(par.M-1) .* (1 - tanE).*(U.^2));
fV_U = taun*(Re.*thSn_prime);
fV_V = -taun;

% Jacobian                
   dE = [par.gamma.*L2 + omega.*L1 + spdiags(fU_U,0,numPar.nx*(numPar.ny-1)+1,numPar.nx*(numPar.ny-1)+1); ...
            spdiags(fV_U,0,numPar.nx*(numPar.ny-1)+1,numPar.nx*(numPar.ny-1)+1)];
            
   dn = [spdiags(fU_V,0,numPar.nx*(numPar.ny-1)+1,numPar.nx*(numPar.ny-1)+1); ...
       par.delta*L2 + omega.*L1 + fV_V.*speye(numPar.nx*(numPar.ny-1)+1,numPar.nx*(numPar.ny-1)+1)];
  
   J = [dE, dn];
   J = sparse(J);        

end
