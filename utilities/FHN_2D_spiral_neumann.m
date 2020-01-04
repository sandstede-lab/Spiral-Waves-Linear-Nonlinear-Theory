function [f,J]  = FHN_2D_spiral_neumann(u,par,numPar,phase_cond)
% Solve for a 2D spiral in the Barkley model on a bounded disk of radius R with Neumann boundary conditions
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

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
ep = par.ep;
delta = par.delta;  % Diffusion coefficient of v-equation
R = par.r2;

[numPar.L1, numPar.L2] = ComputeLinearOperator(par,numPar);
L1 = numPar.L1;
L2 = numPar.L2;
L2 = L2./(R.^2); % Scale the laplacian by the outer radius

U = u(1:nx*ny);
V = u(nx*ny+1:2*nx*ny);
omega = u(end); 

% Select the solution evaluated at R/2
u_phase = U.*phase_cond.pc;
u_phase = u_phase(phase_cond.pc>0); 

% FHN model
fU = (1/ep)*(U.*(1-U).*(U-0.5) - (V+b)/a);
line1 =       L2*U + omega*L1*U + fU;
line2 = delta*L2*V + omega*L1*V + U - V;
line3 = (phase_cond.u_star_th)'*(u_phase - phase_cond.u_star);
f = [line1;line2;line3];

% Jacobian              
if (nargout > 1)
	fU_U = (1/ep)*( (1-U).*(U-0.5) - U.*(U-0.5) + U.*(1-U) );
	fU_V = (1/ep)*(-1/a);
    fV_U = 1;
    fV_V = -1;

    phase_jacob = [sparse(1,nx*(ceil(ny/2) - 1)), phase_cond.u_star_th', sparse(1,nx*floor(ny/2)) ]; % Jacobian of the phase condition
    I = speye(nx*ny,nx*ny);

    dU = [L2 + omega*L1 + spdiags(fU_U,0,nx*ny,nx*ny);
        fV_U.*I;
        phase_jacob];

    dV = [fU_V.*I;
        delta*L2 + omega*L1 + fV_V.*I;
        sparse(1,nx*ny)];

    domega = [L1*U;
        L1*V;
        0];

    J = [dU,dV,domega];      
end
