function [F,J]  = Rossler_2D_spiral_neumann(u,par,numPar,phase_cond)
% Solve for a 2D spiral in the Rossler model on a bounded disk of radius R with Neumann boundary conditions
% Input
%       u: [U; V; W; omega] is initial guess where omega is the free parameter
%       L1: d/d psi, differential opeartor
%       L2: 2D radial Laplacian 
%       par, numPar: system and numerical parameters
%       phase_cond: phase condition for the system
% ** Assumes that the radii for linear operators has been scaled to 1
%
% Output
%       F: function to solve for F(U) = 0
%       J: Jacobian, last column is respect to the free parameter

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;
r = par.r2;

L1 = numPar.L1;
L2 = numPar.L2;

U = u(1:nx*ny);
V = u(nx*ny+1:2*nx*ny);
W = u(2*nx*ny+1:3*nx*ny);
omega = u(end); 

% Select the solution evaluated at R/2
u_phase = U.*phase_cond.pc;
u_phase = u_phase(phase_cond.pc>0); 

% Scale radial laplacian by radius
L2 = L2./(r.^2);

% Reaction terms
fU = -V - W;
fV = U + a*V;
fW = U.*W - c.*W + b;

line1 = delta1.*(L2*U) + omega.*(L1*U) + fU;
line2 = delta2.*(L2*V) + omega.*(L1*V) + fV;
line3 = delta3.*(L2*W) + omega.*(L1*W) + fW;
line4 = (phase_cond.u_star_th)'*(u_phase - phase_cond.u_star);

F = [line1;
    line2;
    line3;
    line4];
   
% Jacobian              
if (nargout > 1)
    
    fU_V = -1;  fU_W = -1;
    fV_U = 1;   fV_V = a;
    fW_U = W;   fW_W = U - c;
    
    I = speye(nx*ny,nx*ny);
    phase_jacob = [sparse(1,numPar.nx*(ceil(numPar.ny/2) - 1)), phase_cond.u_star_th', sparse(1,numPar.nx*floor(numPar.ny/2)) ]; % Only the middle radii
    
    
    dU = [delta1.*L2 + omega.*L1;
        fV_U.*I;
        spdiags(fW_U,0,nx*ny,nx*ny);
        phase_jacob];
    
    dV = [fU_V.*I;
        delta2.*L2 + omega.*L1 + fV_V.*I;
        sparse(nx*ny,nx*ny);
        sparse(1,nx*ny)];
    
    dW = [fU_W.*I;
        sparse(nx*ny,nx*ny);
        delta3.*L2 + omega.*L1 + spdiags(fW_W,0,nx*ny,nx*ny);
        sparse(1,nx*ny)];
    
    domega = [L1*U;
        L1*V;
        L1*W;
        0];
    
        
    J = [dU, dV, dW, domega];
    J = sparse(J);
   
          
end


    
    
    
    
    
    
    