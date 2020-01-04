function [F,J]  = Karma_2D_spiral_neumann(u,par,numPar,phase_cond)
% Solve for a 2D spiral in the Karma model on a bounded disk of radius R with Neumann boundary conditions
% Input
%       u: [U; V; omega] is initial guess where omega is the free parameter
%       L1: d/d psi, differential opeartor
%       L2: 2D radial Laplacian 
%       par, numPar: system and numerical parameters
%       phase_cond: phase condition for the system
% ** Assumes that the radii for linear operators has been scaled to 1
%
% Output
%       F: function to solve for F(U) = 0
%       J: Jacobian, last column is respect to the free parameter

% Rename parameters for convience
tauE = 1/par.tauE;
taun = 1/par.taun;
Re = 1/(1 - exp(-par.Re));
r = par.r2; % Outer radius 

L1 = numPar.L1;
L2 = numPar.L2;

% Scale radial laplacian by radius
L2 = L2./(r.^2);

% Initial guess
U = u(1:numPar.nx*numPar.ny);
V = u(numPar.nx*numPar.ny+1:2*numPar.nx*numPar.ny);
omega = u(end); % Free parameter

% Select the solution evaluated at R/2
u_phase = U.*phase_cond.pc;
u_phase = u_phase(phase_cond.pc>0); 

% Non-linear terms in Karma model
[thSn,thSn_prime] = thetaS(U - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(U - par.Eh);        

fE = tauE .* ( -U + 0.5*(par.Estar - V.^par.M).*(1 - tanE) .* U.^2 );
fn = taun .* ( Re.* thSn - V);    

line1 = par.gamma*(L2*U) + omega*(L1*U) + fE;
line2 = par.delta*(L2*V) + omega*(L1*V) + fn;
line3 = (phase_cond.u_star_th)'*(u_phase - phase_cond.u_star);

F = [line1;line2;line3];
   
% Jacobian              
if (nargout > 1)
    phase_jacob = [sparse(1,numPar.nx*(ceil(numPar.ny/2) - 1)), phase_cond.u_star_th', sparse(1,numPar.nx*floor(numPar.ny/2)) ]; % Jacobian of the phase condition - only acts at R/2
    I = speye(numPar.nx*numPar.ny);
   
    fU_U = tauE*(-1 - 0.5.*(par.Estar - V.^par.M).*(sech(U - par.Eh).^2).* (U.^2) + (par.Estar - V.^par.M).*(1 - tanE).*U);
    fU_V = tauE*(-0.5* par.M.*V.^(par.M-1) .* (1 - tanE).*(U.^2));
    fV_U = taun*(Re.*thSn_prime);
    fV_V = -taun;
   
     
   dE = [par.gamma.*L2 + omega.*L1 + spdiags(fU_U ,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny);
            spdiags(fV_U ,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny) ;
            phase_jacob];
            
   dn = [spdiags(fU_V,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny);
       par.delta*L2 + omega*L1 + fV_V.*I; 
       sparse(1,numPar.nx*numPar.ny)];
   
   domega = [L1*U;
        L1*V;
        0];
   
   J = [dE, dn, domega];
   
   J = sparse(J);
   
          
end


    
    
    
    
    
    
    